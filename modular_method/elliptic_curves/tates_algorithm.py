r"""An implementation of Tate's algorithm for Frey curves

This code allows for the computation of the possible conductors of an
elliptic curve which depends on some integral parameters.

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Joey van Langen (2019-02-27): initial version
- Joey van Langen (2019-10-10): Added termination condition to step 7
- Joey van Langen (2021-02-09): Switched to a queue implementation

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from queue import Queue

from sage.structure.sage_object import SageObject

from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class
from sage.schemes.elliptic_curves.constructor import EllipticCurve

from sage.all import ZZ, QQ, Integer, Infinity
from sage.functions.other import ceil
from sage.misc.functional import is_even, is_odd, norm

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic

from modular_method.padics.pAdic_base import pAdicBase
from modular_method.padics.pAdic_tree import pAdicNode, pAdicTree
from modular_method.padics.pAdic_solver import find_pAdic_roots

from modular_method.diophantine_equations.conditions import TreeCondition, ConditionalValue

from modular_method.elliptic_curves.local_data import FreyCurveLocalData

def tates_algorithm(elliptic_curve, coefficient_ring=None, pAdics=None,
                    base_ring=None, prime=None, initial_values=None,
                    only_calculate=[], precision_cap=Integer(20) , verbose=False):
    r"""Perform Tate's Algorithm on an elliptic curve dependent on
    parameters.

    We consider an elliptic curve $E$ defined by a Weierstrass
    equation of the form ..MATH::

        Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

    where $a_1, a_2, a_3, a_4, a_6$ are (multivariate) polynomials
    over some number field $L$ or a subring thereof. Considering the
    variables of these polynomials to take on undetermined values in
    the ring of integers of a subfield $K$ of $L$ we can see this as
    an elliptic curve over the field $L$.

    To determine the local behavior of $E$ at a finite prime $Q$ of
    $L$ we can follow Tate's algorithm, which in each step makes a
    decision based only on the roots of a polynomial in $a_1, a_2,
    a_3, a_4, a_6$ modulo some power of $Q$. Therefore each step only
    depends on the value of the variables modulo some finite power of
    the prime $P$ in $K$ below $Q$. Using p-adic trees and the
    function :func:`find_pAdic_roots` we can explicitly compute these
    values.

    .. NOTE::

    Even though each step in Tate's algorithm only depends on the
    value of the variables modulo some finite power of $P$ it is not
    guaranteed that the complete algorithm depends only on a finite
    power of $P$. In fact the termination of Tate's algorithm depends
    on the fact that the discriminant of an elliptic curve is bounded.
    This might of course not be the case when the discriminant depends
    on variables.

    In practice the above is only an issue if there is always a case
    in which a new minimal model of lower discriminant exists, or one
    tries to compute all local data and enters the subalgorithm of
    step 7. In all other cases the implementation of this algorithm
    will terminate, but one might enter a curve for which the number
    of steps is very large. It is therefore advised to turn 'verbose'
    on and watch the steps carefully if the computation takes too long
    or prints warnings.
    
    INPUT:
    
    - ``elliptic_curve`` -- The elliptic curve $E$
      
    - ``coefficient_ring`` -- The (multivariate) polynomial ring in
      which $a_1, a_2, a_3, a_4, a_6$ live. By default will be
      initialized as the base ring of $E$.
      
    - ``pAdics`` -- A pAdicBase object or None (default: None). This
      should give the p-adics on $L$ with respect to $Q$. If set to
      None, it will be initialized using the arguments `base_ring` and
      `prime`.
      
    - ``base_ring`` -- A ring or None (default: None). This should be
      a ring whose field of fractions is the field $L$. If set to None
      it will be initialized as the base ring of the argument
      `coefficient_ring`.
      
    - ``prime`` -- The prime $Q$ or None (default None). This may be
      given as a prime ideal or as a generator thereof if it is
      principal. If set to None it will be initialized using the
      p-adics given by the argument `pAdics`. It must be set if that
      argument is set to None.
      
    - ``initial_values`` -- A p-adic tree or None (default:
      None). This should be a tree containing the possible values for
      the variables in the argument `coefficient_ring`. It should be
      given as a pAdicTree that contains the same variables as
      `coefficient_ring`. If set to None, will be initialized as the
      full tree with these variables.
      
    - ``only_calculate`` -- (default: []) When set to a specific value
      or list of values, will ensure that the algorithm only
      calculates a certain quantity or quantities. Possible values are

        - 'conductor' -> Only calculate the conductor exponent

        - 'reduction_type' -> Only calculate the reduction type, i.e.
          good, split multiplicative, non-split multiplicative or
          additive, returned as None, 1, -1 or 0 respectively.

        - 'discriminant' -> Only calculate the valuation of the
          minimal discriminant

        - 'type' -> Only calculate the Kodaira Symbol

        - 'minimal_model' -> Only calculate the minimal model for this
          elliptic curve

        - 'isomorphism' -> Only calculate the change in Weierstrass
          model required to change the given curve into its minimal
          model and the change in Weierstrass model to change the
          minimal model into the given curve.

      By default the function computes all these quantities, but in
      this way they can be selected The function will skip over all
      computations that are not required to determine those quantities
      so this might save time.

    - ``precision_cap`` -- A non-negative integer (default: 20). This
      argument determines the highest precision that will be used for
      the values of the variables of `coefficient_ring`. Note that
      setting this too low might result in inaccurate results, for
      which a warning will appear if this is the case. Setting this
      argument too high might result in endless and slow computations.
    
    - ``verbose`` -- A boolean value or an integer (default: False).
      When set to True or any value larger then zero will print
      comments to stdout about the computations being done whilst
      busy. If set to False or 0 will not print such comments. If set
      to any negative value will also prevent the printing of any
      warnings. Will print more message if set to a higher value.
    
    OUTPUT:
    
    The output will be a FreyCurveLocalData object containing the
    local data of the elliptic curve $E$ for the prime $Q$. If the
    argument `only_calculate` was not the empty list will instead
    return a list consisting of the requested values to be computed in
    the order they were requested in in `only_calculate`.

    If the output depends on the specific values of the variables,
    will return a ConditionalValue with for each case the value as
    mentioned above and the condition that should hold on the
    variables for this value to be attained.
    
    EXAMPLES:
    
    To be made

    """
    _check_elliptic_curve(elliptic_curve)
    coefficient_ring = _init_coefficient_ring(coefficient_ring,
                                              elliptic_curve)
    pAdics = _init_pAdics(pAdics, base_ring, prime, coefficient_ring)
    S = _init_polynomial_ring(coefficient_ring, pAdics)
    variables = _init_variables_tate(S)
    T = _init_initial_values(initial_values, pAdics, variables)
    cases = Queue()
    doneCases = []
    case = _init_case(T, elliptic_curve, S, pAdics, verbose)
    cases.put(case)
    only_calculate = _init_str_list(only_calculate)            

    # Let's work through the queue
    while not cases.empty():
        case = cases.get()
        _tate_step(case, variables, only_calculate, cases,
                   doneCases, verbose, precision_cap)
        cases.task_done()
    # Should consider putting multiple threads on this,
    # but should find a good way of raising exceptions first
        
    return _tate_cleanup(doneCases)

def tates_algorithm_multiple(elliptic_curves, coefficient_rings=None,
                             pAdics=None, base_rings=None, primes=None,
                             initial_values=None, only_calculate=[],
                             input_data=False, single_case=True,
                             precision_cap=Integer(20) , verbose=False):
    r"""Perform Tate's Algorithm on multiple elliptic curves dependent on
    parameters at once.

    Given multiple elliptic curves $E$ as mentioned in
    :func:`tates_algorithm` all depending on the same parameters, will
    do the same computation as :func:`tates_algorithm` for each of
    these at once. Note that the p-adics chosen for each curve may be
    different as long as the p-adics on the parameters is the same for
    each curve.

    One would mainly use this function if for each possible value of
    the parameters it is sufficient to compute the result of Tate's
    algorithm for only one of the given curves.
    
    INPUT:
    
    - ``elliptic_curves`` -- A tuple of elliptic curves $E$
      
    - ``coefficient_rings`` -- A tuple of (multivariate) polynomial
      rings, one for each given curve $E$. Each of these should be the
      ring in which the coefficients of the corresponding curve $E$
      live. By default each of these will be initialized as the base
      ring of the corresponding curve $E$.
      
    - ``pAdics`` -- A tuple of pAdicBase objects, one for each given
      curve $E$. Each entry should be the p-adics to be used for
      Tate's algorithm for the corresponding curve $E$. By default
      this will be initialized using the arguments `base_rings` and
      `primes`.
      
    - ``base_rings`` -- A tuple of rings, one for each given curve
      $E$. The field of fractions of each of these rings should be the
      field over which $E$ is defined when the parameters are
      evaluated. By default each entry will be initialized as the base
      ring of the corresponding entry in `coefficient_rings`.
      
    - ``primes`` -- A tuple of primes, one corresponding to each given
      curve $E$. Each primes should be a finite prime of the fraction
      field of the corresponding entry in `base_ring`. Each prime may
      be given as a prime ideal of the ring of integers or as a
      generator thereof if it is principal. By default it will be
      initialized using the p-adics given by the argument `pAdics`. It
      must be set if that argument is set to None.
      
    - ``initial_values`` -- A p-adic tree or None (default:
      None). This should be a tree containing the possible values for
      the variables in the argument `coefficient_ring`. It should be
      given as a pAdicTree that contains the same variables as all
      rings in `coefficient_rings`. If set to None, will be
      initialized as the full tree with these variables. Note that all
      of the p-adics given by `pAdics` should extend the p-adics of
      this tree.
      
    - ``only_calculate`` -- (default: []) When set to a specific value
      or list of values, will ensure that the algorithm only
      calculates a certain quantity or quantities. Possible values are

        - 'conductor' -> Only calculate the conductor exponent

        - 'reduction_type' -> Only calculate the reduction type, i.e.
          good, split multiplicative, non-split multiplicative or
          additive, returned as None, 1, -1 or 0 respectively.

        - 'discriminant' -> Only calculate the valuation of the
          minimal discriminant

        - 'type' -> Only calculate the Kodaira Symbol

        - 'minimal_model' -> Only calculate the minimal model for this
          elliptic curve

        - 'isomorphism' -> Only calculate the change in Weierstrass
          model required to change the given curve into its minimal
          model and the change in Weierstrass model to change the
          minimal model into the given curve.

      By default the function computes all these quantities, but in
      this way they can be selected. The function will skip over all
      computations that are not required to determine those quantities
      so this might save time.

    - ``input_data`` -- A boolean value (default: False). When set to
      True the output returned will also contain the input data,
      elliptic curve and p-adics, for each case returned by the
      function.

    - ``single_case`` -- A boolean value (default: True). When set to
      True whenever a result is computed for one of the curves, the
      corresponding parameter values are thereafter ignored for the
      computation of all other curves. This ensures the output has a
      unique value for each possible value of the parameters, but
      means this function does not compute Tate's algorithm fully for
      each given curve.

    - ``precision_cap`` -- A non-negative integer (default: 20). This
      argument determines the highest precision that will be used for
      the values of the variables of `coefficient_ring`. Note that
      setting this too low might result in inaccurate results, for
      which a warning will appear if this is the case. Setting this
      argument too high might result in endless and slow computations.
    
    - ``verbose`` -- A boolean value or an integer (default: False).
      When set to True or any value larger then zero will print
      comments to stdout about the computations being done whilst
      busy. If set to False or 0 will not print such comments. If set
      to any negative value will also prevent the printing of any
      warnings. Will print more message if set to a higher value.
    
    OUTPUT:
    
    The output will be a FreyCurveLocalData object containing the
    local data of one of the given curves $E$ at the corresponding
    prime $Q$. If the argument `only_calculate` was not the empty list
    will instead return a list consisting of the requested values to
    be computed in the order they were requested in in
    `only_calculate`.

    If the argument `input_data` was set to `True` the output will be
    a tuple consisting of the input data and the output listed
    above. The input data in this case would be a tuple consisting of
    the specific curve $E$ and the corresponding p-adics to which the
    given output corresponds.

    If the output depends on the specific values of the variables or
    the given elliptic curve, will return a ConditionalValue with for
    each case the value as mentioned above and the condition that
    should hold on the variables for this value to be attained. Note
    that if `single_case` was set to `False` this ConditionalValue
    might have overlapping conditions as each given elliptic curve $E$
    may have a different outcome.

    """
    n = len(elliptic_curves)
    if coefficient_rings is None:
        coefficient_rings = [None for i in range(n)]
    else:
        coefficient_rings = list(coefficient_rings)
    if pAdics is None:
        pAdics = [None for i in range(n)]
    else:
        pAdics = list(pAdics)
    if base_rings is None:
        base_rings = tuple(None for i in range(n))
    if primes is None:
        primes = tuple(None for i in range(n))
    cases = Queue()
    doneCases = []
    for i in range(n):
        _check_elliptic_curve(elliptic_curves[i])
        coefficient_rings[i] = _init_coefficient_ring(coefficient_rings[i],
                                                      elliptic_curves[i])
        pAdics[i] = _init_pAdics(pAdics[i], base_rings[i], primes[i],
                                 coefficient_rings[i])
        S = _init_polynomial_ring(coefficient_rings[i], pAdics[i])
        variables = _init_variables_tate(S)
        initial_values = _init_initial_values(initial_values,
                                              pAdics[i], variables)
        case = _init_case(initial_values, elliptic_curves[i], S,
                          pAdics[i], verbose)
        case['disjoint'] = Integer(0) 
        if input_data:
            case['Estart'] = case['E']
        cases.put(case)
    only_calculate = _init_str_list(only_calculate)            

    # Let's work through the queue
    while not cases.empty():
        case = cases.get()
        i = case['disjoint']
        do_case = True
        while (single_case and i < len(doneCases)):
            case['T'].cut(doneCases[i][Integer(1) ]._root)
            case['disjoint'] += Integer(1) 
            i = case['disjoint']
            if case['T'].is_empty(): # Unnecessary, so skip
                if verbose:
                    if verbose > Integer(1) :
                        print(f"Skipped the case {case['number']} as it is " +
                              "not wanted for any parameter values anymore")
                    else:
                        print(f"Skipped a case as it is not wanted " +
                              "for any parameter values anymore")
                do_case = False
                break
        if do_case:
            _tate_step(case, variables, only_calculate, cases, doneCases,
                       verbose, precision_cap, input_data=input_data)
        cases.task_done()
    # Should consider putting multiple threads on this,
    # but should find a good way of raising exceptions first
        
    return _tate_cleanup(doneCases)

def _tate_step(case, variables, only_calculate, cases, doneCases,
               verbose, precision_cap, input_data=False):
    r"""Perform a single step of Tate's algorithm on a given case."""
    if 'next' in case:
        if case['next'] == 'Step1':
            if verbose > Integer(0) :
                printstr = "Performing step 1"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step1(case, cases, variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step2T':         
            if verbose > Integer(0) :
                printstr = "Performing the transformation for step 2"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step2_t(case, cases,
                          verbose=(verbose and verbose-Integer(1) ))
        elif case['next'] == 'Step2':
            if verbose > Integer(0) :
                printstr = "Performing step 2"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step2(case, cases,
                        variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step3':
            if verbose > Integer(0) :
                printstr = "Performing step 3"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step3(case, cases,
                        variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step4':
            if verbose > Integer(0) :
                printstr = "Performing step 4"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step4(case, cases,
                        variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step5':
            if verbose > Integer(0) :
                printstr = "Performing step 5"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step5(case, cases,
                        variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step6T':
            if verbose > Integer(0) :
                printstr = "Performing the transformation for step 6"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step6_t(case, cases,
                          verbose=(verbose and verbose-Integer(1) ))
        elif case['next'] == 'Step6':
            if verbose > Integer(0) :
                printstr = "Performing step 6"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step6(case, cases,
                        variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step7':
            if verbose > Integer(0) :
                printstr = "Performing step 7"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step7(case, cases, only_calculate,
                        variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'].startswith('Step7('):
            if case['next'].endswith('T'):
                n = int(case['next'][Integer(6) :-Integer(2) ])
                if verbose > Integer(0) :
                    printstr = f"Performing the transformation for step 7({n})"
                    if verbose > Integer(1) :
                        printstr += (" for case " +
                                     case['number'])
                    print(printstr)
                _tate_step7sub_t(n, case, cases,
                                 verbose=(verbose and verbose-Integer(1) ))
            else:
                n = int(case['next'][Integer(6) :-Integer(1) ])
                if verbose > Integer(0) :
                    printstr = f"Performing step 7({n})"
                    if verbose > Integer(1) :
                        printstr += (" for case " +
                                     case['number'])
                    print(printstr)
                _tate_step7sub(n, case, cases, only_calculate,
                               variables=variables,
                               verbose=(verbose and verbose-Integer(1) ),
                               precision_cap=precision_cap)
        elif case['next'] == 'Step8T':
            if verbose > Integer(0) :
                printstr = "Performing the transformation for step 8"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step8_t(case, cases,
                          verbose=(verbose and verbose-Integer(1) ))
        elif case['next'] == 'Step8':
            if verbose > Integer(0) :
                printstr = "Performing step 8"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step8(case, cases, variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step9T':
            if verbose > Integer(0) :
                printstr = "Performing the transformation for step 9"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step9_t(case, cases,
                          verbose=(verbose and verbose-Integer(1) ))
        elif case['next'] == 'Step9':
            if verbose > Integer(0) :
                printstr = "Performing step 9"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step9(case, cases, variables=variables,
                        verbose=(verbose and verbose-Integer(1) ),
                        precision_cap=precision_cap)
        elif case['next'] == 'Step10':
            if verbose > Integer(0) :
                printstr = "Performing step 10"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step10(case, cases, variables=variables,
                         verbose=(verbose and verbose-Integer(1) ),
                         precision_cap=precision_cap)
        elif case['next'] == 'Step11':
            if verbose > Integer(0) :
                printstr = "Performing step 11"
                if verbose > Integer(1) :
                    printstr += (" for case " +
                                 case['number'])
                print(printstr)
            _tate_step11(case, cases,
                         verbose=(verbose and verbose-Integer(1) ))
        else:
            print(f"Skipping case for unknown step: {case['next']}")
    elif _should_calculate_vDelta(case, only_calculate):
        if verbose > Integer(0) :
            printstr = "Calculating valuation of discriminant"
            if verbose > Integer(1) :
                printstr += (" for case " +
                             case['number'])
            print(printstr)
        _tate_calculate_vDelta(case, cases, variables=variables,
                               verbose=(verbose and verbose-Integer(1) ),
                               precision_cap=precision_cap)
    elif _should_calculate_m(case, only_calculate):
        if verbose > Integer(0) :
            printstr = "Calculating number of components of the special fiber"
            if verbose > Integer(1) :
                printstr += (" for case " +
                             case['number'])
            print(printstr)
        _tate_calculate_m(case, cases)
    elif _should_calculate_f(case, only_calculate):
        if verbose > Integer(0) :
            printstr = "Calculating exponent of the conductor"
            if verbose > Integer(1) :
                printstr += (" for case " +
                             case['number'])
            print(printstr)
        _tate_calculate_f(case, cases)
    elif _should_calculate_n(case, only_calculate):
        if verbose > Integer(0) :
            printstr = "Calculating the index of the Kodaira symbol"
            if verbose > Integer(1) :
                printstr += (" for case " +
                             case['number'])
            print(printstr)
        _tate_calculate_n(case, cases)
    elif _should_calculate_split(case, only_calculate):
        if verbose > Integer(0) :
            printstr = "Calculating whether multiplicative reduction is split"
            if verbose > Integer(1) :
                printstr += (" for case " +
                             case['number'])
            print(printstr)
        _tate_calculate_split(case, cases, variables=variables,
                              verbose=(verbose and verbose-Integer(1) ),
                              precision_cap=precision_cap)
    elif _should_calculate_c(case, only_calculate):
        if verbose > Integer(0) :
            printstr = "Calculating the order of the group of components"
            if verbose > Integer(1) :
                printstr += (" for case " +
                             case['number'])
            print(printstr)
        _tate_calculate_c(case, cases, variables=variables,
                          verbose=(verbose and verbose-Integer(1) ),
                          precision_cap=precision_cap)
    else:
        if verbose > Integer(0) :
            if verbose > Integer(1) :
                print("Finishing case " +
                      case['number'])
            else:
                print("Finishing a case")
        _tate_finish(case, only_calculate, result=doneCases,
                     variables=variables, input_data=input_data)

# An easy way to keep track of cases when printing verbose
_case_number = Integer(0) 

def _get_case_number():
    global _case_number
    _case_number += Integer(1) 
    return str(_case_number)

def _init_case_number():
    global _case_number
    _case_number = Integer(0) 

def _least_power(poly_list, pAdics):
    r"""Get the smallest valuation among coefficients of polynomials

    INPUT:

    - ``poly_list`` -- A list of polynomials

    - ``pAdics`` -- The p-adics to be used for the valuation

    OUTPUT:

    The smallest valuation among the coefficients of the given
    polynomials, or Infinity if there were no non-zero polynomials.

    """
    result = Infinity
    for poly in poly_list:
        result = min([result]+[pAdics.valuation(c)
                               for c in poly.coefficients()])
    return result

def _determine_level(poly_list, pAdics, T, max_value):
    r"""Determine the level of nodes sufficient for polynomials to have a
    fixed value modulo some prime power.

    Let $L / K$ be an extension of number fields and $Q$ be a finite
    prime of $L$ lying above the prime $P$ of $K$. For a list of
    (multivariate) polynomials $f_1, \ldots, f_n$ over $L$ one can
    wonder up to what power $r$ of the prime $P$ we have to determine
    the value of the variables in the ring of integers of $K$ to fix
    the value of all polynomials modulo $Q^s$. Given $s$ this function
    computes $r$.

    INPUT:

    - ``poly_list`` -- The list of polynomials $f_1, \ldots, f_n$
    
    - ``pAdics`` -- The p-adics determined by $L$ and $Q$

    - ``T`` -- The root of the p-adic tree containing the possible
      values for the variables

    - ``max_value`` -- The integer $s$

    OUTPUT:

    The smallest possible integer $r$ or $0$ if it does not exist or
    is smaller than $0$.

    """
    level = max_value - _least_power(poly_list, pAdics)
    if level > -Infinity:
        level = ceil(level/pAdics.extension_multiplicity(T.pAdics()))
    if level < Integer(0) :
        level = Integer(0) 
    return level

def _get_cases_invariant(poly, pAdics, T, name, general_case, result,
                         variables=None, verbose=False, precision=Integer(20) ,
                         precision_cap=Integer(20) , **kwds):
    r"""Determine the different valuations a polynomial can have.

    For each possible valuation of the polynomial `poly`, puts a copy
    of `general_case` in `result`. Each copy has two additional
    entries: An entry with key `name` and value a possible valuation
    of the polynomial and an entry with key 'T' and value the root of
    a p-adic tree containing those values in the given p-adic tree for
    which the corresponding valuation of the polynomial is attained.

    INPUT:

    - ``poly`` -- A (multivariate) polynomial

    - ``pAdics`` -- The p-adics to be used

    - ``T`` -- The root of the p-adic tree containing the possible
      values for the variables

    - ``name`` -- A string which is the name to be assigned to the
      valuation of the polynomial

    - ``general_case`` -- A dictionary that is a template for the
      cases to be returned.

    - ``result`` -- A queue in which to put the different cases of
      this computation

    - ``variables`` -- A list of the variables of the polynomial or
      None (default) for it to be initialized from the polynomial.

    - ``verbose`` -- Verbosity argument

    - ``precision`` -- The precision to be used when computing the
      valuation

    - ``precision_cap`` -- The largest precision to be used on the
      variables

    TESTS:

    In case the actual valuation is above the given precision, this
    function should change the precision to find the valuation, only
    stopping for the given precision_cap::

        sage: from queue import Queue
        sage: from modular_method.elliptic_curves.tates_algorithm import _get_cases_invariant
        sage: from modular_method.padics.pAdic_base import pAdicBase
        sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
        sage: R.<x> = QQ[]
        sage: T = CongruenceCondition(x, 2).pAdic_tree(pAdics=pAdicBase(QQ, 2))
        sage: N = 2^18*x^8 + 2^24
        sage: queue = Queue()
        sage: _get_cases_invariant(N, T.pAdics(), T.root(), 'Nval', {}, queue, variables=[x])
        sage: result = queue.get_nowait()
        sage: queue.empty()
        True
        sage: result['Nval']
        24
        sage: result['T'] == T.root()
        True

    """
    Tlist, v_min = find_pAdic_roots(poly, pAdics=pAdics,
                                    variables=variables, value_tree=T,
                                    precision=precision,
                                    verbose=(verbose-Integer(1)  if verbose > Integer(0) 
                                             else verbose),
                                    give_list=True,
                                    precision_cap=precision_cap)
    while (v_min + len(Tlist) > precision and # Happens when precision was
           not Tlist[-Integer(1) ].is_empty()):         # reached, should refine search
        precision = precision + Integer(10) 
        Tlist_, v_min_ = find_pAdic_roots(poly, pAdics=pAdics,
                                          variables=variables,
                                          value_tree=Tlist.pop(),
                                          precision=precision,
                                          verbose=(verbose-Integer(1)  if
                                                   verbose > Integer(0)  else verbose),
                                          give_list=True,
                                          precision_cap=precision_cap)
        i = v_min - v_min_ + len(Tlist)
        Tlist.extend(Tlist_[i:])
    for i in range(len(Tlist)):
        if not Tlist[i].is_empty():
            case = general_case.copy()
            case[name] = v_min + i
            case['T'] = Tlist[i]   
            if verbose > Integer(0) :
                case['number'] = _get_case_number()
                print("Added case " + str(case['number']) + ": valuation of")
                print(poly)
                print("equal to " + str(v_min + i) + " for ")
                n = Tlist[i].minimum_full_level()
                resultstr = ""
                for N in Tlist[i].children_at_level(n):
                    resultstr += str(N.representative()) + " "
                resultstr += ("modulo " + str(Tlist[i].pAdics().prime()) +
                           "^" + str(n))
                print(resultstr)
            result.put(case)

def _get_two_cases_invariant(poly, pAdics, T, bndry, case_big,
                             case_small, result, variables=None,
                             verbose=False, precision_cap=Integer(20) , **kwds):
    r"""Distinguish whether a polynomial has at least a certain valuation
    or not

    In case there is possible values for the variables in which the
    valuation of the polynomial is at least `bndry` then the case
    `case_big` is put in `result` with an added entry with key 'T' and
    value the root of a p-adic tree containing all those values for
    which this is the case.

    If there is possible values for the variables in which the
    valuation of the polynomial is less then `bndry` then the case
    `case_small` is put in `result` with and added entry with key 'T'
    and value the root of a p-adic tree containing all those values
    for which this is the case.

    INPUT:

    - ``poly`` -- The polynomial

    - ``pAdics`` -- The p-adics to be used

    - ``T`` -- The root of a p-adic tree containing the possible
      values for the variables

    - ``bndry`` -- An integer giving the valuation at which we
      distinguish

    - ``case_big`` -- A dictionary that is a template for the case the
      valuation is at least `bndry`

    - ``case_small`` -- A dictionary that is a template for the case
      the valuation is less then `bndry`

    - ``result`` -- A queue in which to put the different cases of
      this computation

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from the polynomial
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    Tbig, Tsmall = find_pAdic_roots(poly, pAdics=pAdics, variables=variables,
                                    value_tree=T, precision=bndry,
                                    verbose=(verbose-Integer(1)  if verbose > Integer(0) 
                                             else verbose),
                                    precision_cap=precision_cap)
    if not Tbig.is_empty():
        case_big['T'] = Tbig
        if verbose > Integer(0) :
            case_big['number'] = _get_case_number()
            print("Added case " + str(case_big['number']) +
                  ": valuation of")
            print(poly)
            print("at least " + str(bndry) + " for ")
            n = Tbig.minimum_full_level()
            resultstr = ""
            for N in Tbig.children_at_level(n):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(Tbig.pAdics().prime()) +
                       "^" + str(n))
            print(resultstr)
        result.put(case_big)
    if not Tsmall.is_empty():
        case_small['T'] = Tsmall
        if verbose > Integer(0) :
            case_small['number'] = _get_case_number()
            print("Added case " + str(case_small['number']) +
                  ": valuation of")
            print(poly)
            print("less than " + str(bndry) + " for ")
            n = Tsmall.minimum_full_level()
            resultstr = ""
            for N in Tsmall.children_at_level(n):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(Tsmall.pAdics().prime()) +
                          "^" + str(n))
            print(resultstr)
        result.put(case_small)
        
def _get_number_of_roots_cases(poly_list, pAdics, T, name,
                               general_case, result, variables=None,
                               verbose=False, **kwds):
    r"""Determine the possible roots a polynomial with as coefficients
    multivariate polynomials can have.

    Let $L / K$ be an extension of number fields and let $Q$ be a
    finite prime of $L$ above the prime $P$ of $K$. For (multivariate)
    polynomials $f_0, \ldots, f_n$ over $L$ we can construct the
    polynomial .. MATH::
    
        F(X) = f_0*X^n + ... + f_{n-1}*X + f_n

    If the variables of the polynomials $f_0, \ldots, f_n$ take on
    values in the ring of integers of $K$ such that the valuations of
    $f_0, \ldots, f_n$ are non-negative, the above defines a
    polynomial over the residue field of $Q$. This function determines
    how many roots this polynomial has over that residue field and for
    which value that number of roots occurs.

    For each possible number of roots of the polynomial $F$ over the
    residue field of $Q$ a copy of the case `general_case` is put into
    `result` with as additional entries:

    - An entry with as key the string given by `name` and as value a
      number of roots the polynomial $F$ has over the residue field of
      $Q$.

    - An entry with as key 'T' and as value the root of a p-adic tree
      containing all those values for the variables in the given tree
      for which this number of roots for $F$ is attained.

    INPUT:

    - ``poly_list`` -- A list of the polynomials $f_0, \ldots, f_n$ in
      that order

    - ``pAdics`` -- The p-adics defined by $L$ and $Q$

    - ``T`` -- The root of a p-adic tree giving the possible values
      for the variables. These should be such that the valuation of
      $f_0, \ldots, f_n$ should be non-negative.
    
    - ``name`` -- A string that is the name that should be assigned to
      the number of roots of the polynomial $F$

    - ``general_case`` -- A dictionary that is a template for the
      cases to be returned

    - ``result`` -- A queue in which to put the different cases of
      this computation

    - ``variables`` -- A list of the variables of the polynomials or
      None if it should be determined from the polynomials

    - ``verbose`` -- Verbosity argument

    """
    Tdict = {}
    n = _determine_level(poly_list, pAdics, T, Integer(1) )
    F = pAdics.residue_field()
    S = F['x']; (x,) = S._first_ngens(1)
    child_list = T.children_at_level(n)
    if verbose > Integer(0) :
        print("Checking irreducibility in %d cases"%len(child_list))
    for child in child_list:
        coeff_list = [F(poly(child.representative())) for poly in poly_list]
        f = Integer(0) 
        k = len(coeff_list)-Integer(1) 
        for i in range(k+Integer(1) ):
            f += coeff_list[i] * x**(k-i)
        m = sum([r[Integer(1) ] for r in f.roots()])
        if not (m in Tdict):
            Tdict[m] = pAdicNode(pAdics=T.pAdics(), width=T.width)
        Tdict[m].merge(child, from_root=True)
    for (m,Tm) in Tdict.items():
        case = general_case.copy()
        case[name] = m
        case['T'] = Tm
        if verbose > Integer(0) :
            case['number'] = _get_case_number()
            print("Added case " + str(case['number']) + ":")
            resultstr = ""
            for i in range(len(poly_list)):
                resultstr += "(" + str(poly_list[i]) + ")"
                if i < len(poly_list):
                    resultstr += (" * T^" +
                                  str(len(poly_list)-i) +
                                  " + ")
            print(resultstr)
            print("has exactly " + str(m) + " roots for ")
            nT = Tm.minimum_full_level()
            resultstr = ""
            for N in Tm.children_at_level(nT):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(Tm.pAdics().prime()) +
                       "^" + str(nT))
            print(resultstr)
        result.put(case)
    
def _tate_step1(case, result, **kwds):
    r"""Perform step 1 of Tate's algorithm.

    Stop if the valuation of the discriminant of the elliptic curve is
    less than 1. Continue to step 2 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $I_0$

    - valuation of the minimal discriminant: 0
    
    - number of components: 1

    - exponent of the conductor: 0

    - order of the group of components: 1

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    case_big = case.copy()
    case_big['next'] = 'Step2T'
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(vDelta=Integer(0) , m=Integer(1) , f=Integer(0) , c=Integer(1) , KS="I0"))
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    _get_two_cases_invariant(S(E.discriminant()), pAdics, T, Integer(1) ,
                             case_big, case_small, result, **kwds)
    
def _tate_step2_t(case, result, verbose=False, **kwds):
    r"""Perform the transformation necessary before step 2 of Tate's
    algorithm.

    Move the singular point to (0,0), i.e. transform the curve such
    that it has a Weierstrass equation of the form .. MATH::

        Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

    for which $a_3$, $a_4$ and $a_6$ have valuation at least 1.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step
    
    - ``verbose`` -- Verbosity argument

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    urst = case['urst']
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    char = pAdics.characteristic()
    replaceCases = dict()
    
    # The parameter s is used to compute
    # - square roots in characteristic 2
    # - cube roots in characteristic 3
    # In both cases these roots are unique and can be computed
    # by taking an element to the power s
    s = Integer(0) 
    if char == Integer(2) :
        if R == ZZ and P in ZZ.ideal_monoid():
            s = ((ZZ.quotient(P.gen() - Integer(1) )(Integer(2) ))**(-Integer(1) )).lift()
        else:
            s = ((ZZ.quotient(ZZ(norm(P) - Integer(1) ))(Integer(2) ))**(-Integer(1) )).lift()
    if char == Integer(3) :
        if R == ZZ and P in ZZ.ideal_monoid():
            s = ((ZZ.quotient(P.gen() - Integer(1) )(Integer(3) ))**(-Integer(1) )).lift() 
        else:
            s = ((ZZ.quotient(ZZ(norm(P) - Integer(1) ))(Integer(3) ))**(-Integer(1) )).lift()
    if s == Integer(0) :
        s = Integer(1) 

    # Determining the necessary transformations
    level = _determine_level([S(E.a1()), S(E.a2()), S(E.a3()), S(E.a4()),
                              S(E.a6())], pAdics, T, Integer(1) )
    if verbose > Integer(0) :
        print("Determining singular point for " +
               str(T.count_children_at_level(level)) + " cases")
    for node in T.children_at_level(level):
        # Determining the singular point x,y
        a1 = F(S(E.a1())(node.representative()))
        a2 = F(S(E.a2())(node.representative()))
        a3 = F(S(E.a3())(node.representative()))
        a4 = F(S(E.a4())(node.representative()))
        a6 = F(S(E.a6())(node.representative()))    
        if char == Integer(2) :
            if a1 == Integer(0) :
                x = a4**s
                y = (a6 + a2*a4)**s
            else:
                x = a3 / a1
                y = (x**Integer(2)  + a4) / a1
        else:
            # Coordinate transformation to end up in the case
            # where a1 = a3 = 0, hence y = 0
            a22 = a2 + a1**Integer(2)  / Integer(4) 
            a42 = a4 + a1*a3 / Integer(2) 
            a62 = a6 + a3**Integer(2)  / Integer(4) 
            y = Integer(0) 
            if char == Integer(3) :
                if a22 == Integer(0) :
                    x = (-a62)**s
                else:
                    x = -a42 / a22
            else:
                # Coordinate transformation to end up in the case where
                # also a2 = 0
                a43 = a42 - a22**Integer(2)  / Integer(3) 
                a63 = a62 - a22*a42 / Integer(3)  + Integer(2) *a22**Integer(3)  / Integer(27)                 
                if a43 == Integer(0) :
                    x = Integer(0) 
                else:
                    x = -Integer(3) *a63 / (Integer(2) *a43)                
                # Transforming back
                x = x - a22 / Integer(3)                 
            # Transforming back
            y = y - a1*x / Integer(2)  - a3 / Integer(2) 
        singularPoint = tuple([x,y])
        if singularPoint in replaceCases:
            replaceCases[singularPoint].merge(node, from_root=True)
        else:
            Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
            Tn.merge(node, from_root=True)
            replaceCases[singularPoint] = Tn

    # Doing the actual transformations
    for (point,Tn) in replaceCases.items():
        xn = F.lift(point[Integer(0) ])
        yn = F.lift(point[Integer(1) ])
        En = E.rst_transform(xn,Integer(0) ,yn)
        newCase = case.copy()
        newCase.update(dict(next='Step2', T=Tn, E=En,
                            urst=_urst_combine(urst, (Integer(1) , xn, Integer(0) , yn))))
        if verbose > Integer(0) :
            newCase['number'] = _get_case_number()
            print("Added case " + str(newCase['number']) + ":")
            print("Transform curve into")
            print(En)
            print("for")
            nT = Tn.minimum_full_level()
            resultstr = ""
            for N in Tn.children_at_level(nT):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(Tn.pAdics().prime()) +
                          "^" + str(nT))
            print(resultstr)
        result.put(newCase)
    
def _tate_step2(case, result, **kwds):
    r"""Perform step 2 of Tate's algorithm.

    Stop if the valuation of the invariant $b_2$ of the elliptic curve
    is less than 1. Continue to step 3 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $I_n$

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    case_big = case.copy()
    case_big.update(dict(next='Step3'))
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(f=Integer(1) , KS="In"))
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    _get_two_cases_invariant(S(E.b2()), pAdics, T, Integer(1) , case_big,
                             case_small, result, **kwds)

def _tate_step3(case, result, **kwds):
    r"""Perform step 3 of Tate's algorithm.

    Stop if the valuation of the invariant $a_6$ of the elliptic curve
    is less than 2. Continue to step 4 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $II$
    
    - number of components: 1

    - order of the group of components: 1

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    case_big = case.copy()
    case_big.update(dict(next='Step4'))
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(KS="II", m=Integer(1) , c=Integer(1) ))
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    char = pAdics.characteristic()
    if char != Integer(2)  and char != Integer(3) :
        case_small['vDelta'] = Integer(2) 
        case_small['f'] = Integer(2) 
    _get_two_cases_invariant(S(E.a6()), pAdics, T, Integer(2) , case_big,
                             case_small, result, **kwds)
    
def _tate_step4(case, result, **kwds):
    r"""Perform step 4 of Tate's algorithm.

    Stop if the valuation of the invariant $b_8$ of the elliptic curve
    is less than 3. Continue to step 5 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $III$,
    
    - number of components: 2

    - order of the group of components: 2

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    case_big = case.copy()
    case_big['next'] = 'Step5'
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(KS="III", m=Integer(2) , c=Integer(2) ))
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    char = pAdics.characteristic()
    if char != Integer(2) :
        case_small['vDelta'] = Integer(3) 
        case_small['f'] = Integer(2) 
    _get_two_cases_invariant(S(E.b8()), pAdics, T, Integer(3) , case_big,
                             case_small, result, **kwds)

def _tate_step5(case, result, **kwds):
    r"""Perform step 5 of Tate's algorithm.

    Stop if the valuation of the invariant $b_6$ of the elliptic curve
    is less than 3. Continue to step 6 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $IV$
    
    - number of components: 3

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    case_big = case.copy()
    case_big['next'] = 'Step6T'
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(KS="IV", m=Integer(3) ))
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    char = pAdics.characteristic()
    if char != Integer(3) :
        case_small['vDelta'] = Integer(4) 
        case_small['f'] = Integer(2) 
    _get_two_cases_invariant(S(E.b6()), pAdics, T, Integer(3) , case_big,
                             case_small, result, **kwds)

def _tate_step6_t(case, result, verbose=False):
    r"""Perform the transformation necessary before step 6 of Tate's
    algorithm.

    Transform the curve such that it has a Weierstrass equation of the
    form .. MATH::

        Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

    for which $a_1$ and $a_2$ having valuation at least 1, $a_3$ and
    $a_4$ having valuation at least 2 and $a_6$ having valuation at
    least 3

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``verbose`` -- Verbosity argument

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    urst = case['urst']
    pi = pAdics.uniformizer()
    char = pAdics.characteristic()
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()

    # Determining the diferent transformations needed
    changeDict = dict()   
    if char != Integer(2) :
        a1 = S(E.a1())
        a31 = S(E.a3()/pi)     
        half = (F(Integer(2) )**(-Integer(1) ))
        level = _determine_level([a1, a31], pAdics, T, Integer(1) )
        if verbose > Integer(0) :
            print("Determining necessary transformation for " +
                   str(T.count_children_at_level(level)) + " cases")
        for node in T.children_at_level(level):
            alpha = -half * F(a1(node.representative()))
            a31_ev = F(a31(node.representative()))
            beta = -half * a31_ev
            alphaBetaPair = tuple([alpha, beta])
            if alphaBetaPair in changeDict:
                changeDict[alphaBetaPair].merge(node, from_root=True)
            else:
                T0 = pAdicNode(pAdics=T.pAdics(), width=T.width)
                T0.merge(node, from_root=True)
                changeDict[alphaBetaPair] = T0
    else:
        if R == ZZ and P in ZZ.ideal_monoid():
            sqrtPower = ((ZZ.quotient(P.gen() - Integer(1) )(Integer(2) ))**(-Integer(1) )).lift()
        else:
            sqrtPower = ((ZZ.quotient(ZZ(norm(P) - Integer(1) ))(Integer(2) ))**(-Integer(1) )).lift()
        if sqrtPower == Integer(0) :
            sqrtPower = Integer(1)             
        a2 = S(E.a2())
        a62 = S(E.a6()/(pi**Integer(2) ))        
        level = _determine_level([a2, a62], pAdics, T, Integer(1) )
        if verbose > Integer(0) :
            print("Determining necessary transformation for " +
                   str(T.count_children_at_level(level)) + " cases")
        for node in T.children_at_level(level):
            alphaSqrd = -F(a2(node.representative()))
            betaSqrd = F(-a62(node.representative()))
            alpha = alphaSqrd**sqrtPower
            beta = betaSqrd**sqrtPower
            alphaBetaPair = tuple([alpha, beta])
            if alphaBetaPair in changeDict:
                changeDict[alphaBetaPair].merge(node, from_root=True)
            else:
                Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
                Tn.merge(node, from_root=True)
                changeDict[alphaBetaPair] = Tn

    # Performing the necessary transformations
    for (alphaBetaPair, Tn) in changeDict.items():
        s = F.lift(alphaBetaPair[Integer(0) ])
        t = F.lift(alphaBetaPair[Integer(1) ])*pi
        En = E.rst_transform(Integer(0) , s, t)
        newCase = case.copy()
        newCase.update(dict(next='Step6', T=Tn, E=En,
                            urst=_urst_combine(urst, (Integer(1) , Integer(0) , s, t))))
        if verbose > Integer(0) :
            newCase['number'] = _get_case_number()
            print("Added case " + str(newCase['number']) + ":")
            print("Transform curve into")
            print(En)
            print("for")
            nT = Tn.minimum_full_level()
            resultstr = ""
            for N in Tn.children_at_level(nT):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(Tn.pAdics().prime()) +
                          "^" + str(nT))
            print(resultstr)
        result.put(newCase)

def _tate_step6(case, result, **kwds):
    r"""Perform step 6 of Tate's algorithm.

    Stop if the polynomial .. MATH::
    
        X^3 + \frac{a_2}{\pi} X^2 + \frac{a_4}{\pi^2} X +
        \frac{a_6}{\pi^3}

    has three distinct roots in the residue field of the p-adics. Here
    $a_2$, $a_4$ and $a_6$ are the respective invariants of the
    elliptic curve and $\pi$ is the uniformizer of the
    p-adics. Equivalently stop if the valuation of .. MATH::

        -4 a_2^3 a_6 + a_2^2 a_4^2 - 4 a_4^3 - 27 a_6^2 +
        18 a_2 a_4 a_6

    is less than 7. Otherwise continue to step 7.

    In case we stop we have:
    
    - Kodaira symbol: $I_0^*$
    
    - number of components: 5

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    case_big = case.copy()
    case_big['next'] = 'Step7'
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(KS="I0*", m=Integer(5) ))
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    char = pAdics.characteristic()
    if char != Integer(2) :
        case_small['vDelta'] = Integer(6) 
        case_small['f'] = Integer(2) 
    D = S(-Integer(4) *E.a2()**Integer(3) *E.a6() + E.a2()**Integer(2) *E.a4()**Integer(2)  - Integer(4) *E.a4()**Integer(3)  - Integer(27) *E.a6()**Integer(2)  +
          Integer(18) *E.a2()*E.a4()*E.a6())
    _get_two_cases_invariant(D, pAdics, T, Integer(7) , case_big, case_small,
                             result, **kwds)
        
def _tate_step7(case, result, restrictions, **kwds):
    r"""Perform step 7 of Tate's algorithm.

    Stop if the polynomial .. MATH::
    
        X^3 + \frac{a_2}{\pi} X^2 + \frac{a_4}{\pi^2} X +
        \frac{a_6}{\pi^3}

    has one single and one double root in the residue field of the
    p-adics. Here $a_2$, $a_4$ and $a_6$ are the respective invariants
    of the elliptic curve and $\pi$ is the uniformizer of the
    p-adics. Equivalently stop if the valuation of .. MATH::

        3 a_4 - a_2^2

    is less than 3. Otherwise continue to step 8.

    In case we stop we have:
    
    - Kodaira symbol: $I_n^*$

    We also need to perform the subalgorithm for step 7 if the
    characteristic is 2 or we want the order of the group of
    components.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``restrictions`` -- A list of values to which the computation
      should be limited. See the argument `only_calculate` in
      :func:`tates_algorithm` for more information.

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    case_big = case.copy()
    case_big['next'] = 'Step8T'
    case_small = case.copy()
    if (pAdics.characteristic() == Integer(2)  or 'type' in restrictions or
        _should_calculate_c(case, restrictions)):
        case_small['next'] = 'Step7(1)T'
    else:
        del case_small['next']
        case_small.update(dict(KS="In*", f=Integer(2) ))
    _get_two_cases_invariant(S(Integer(3) *E.a4() - (E.a2())**Integer(2) ), pAdics, T, Integer(3) ,
                             case_big, case_small, result,
                             **kwds)
            
def _tate_step7sub_t(n, case, result, verbose=False):
    r"""Perform the transformation necessary before a step in the
    subalgorithm of step 7 of Tate's algorithm.

    If this is the $n$-th step in the subalgorithm transform the curve
    such that it has a Weierstrass equation of the form .. MATH::

        Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

    for which $a_2$ has valuation precisely 1, $a_3$ has valuation at
    least $(n+4)/2$ if $n$ is even, $a_4$ has valuation at least
    $(n+5)/2$ if $n$ is odd, and $a_6$ has valuation at least $n + 3$.

    INPUT:

    - ``n`` -- The number of the step in the subalgorithm that will be
      performed after this transformation.

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``verbose`` -- Verbosity argument

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    urst = case['urst']
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    pi = pAdics.uniformizer()
    char = pAdics.characteristic()

    # For characteristic 2 compute the integer s
    # such that a^s is the square root of a
    if char == Integer(2) :
        if R == ZZ and P in ZZ.ideal_monoid():
            sqrtPower = ((ZZ.quotient(P.gen() - Integer(1) )(Integer(2) ))**(-Integer(1) )).lift()
        else:
            sqrtPower = ((ZZ.quotient(ZZ(norm(P) - Integer(1) ))(Integer(2) ))**(-Integer(1) )).lift()
        if sqrtPower == Integer(0) :
            sqrtPower = Integer(1) 

    # Determine the level of nodes responsible for the transformations
    changeDict = dict()
    if n == Integer(1) :
        a42 = S(E.a4()/(pi**Integer(2) ))
        if char == Integer(2) :
            level = _determine_level([a42], pAdics, T, Integer(1) )
        else:
            a21 = S(E.a2()/pi)
            if char == Integer(3) :
                level = _determine_level([a21, a42], pAdics, T, Integer(1) )
            else:
                a63 = S(E.a6()/(pi**Integer(3) ))
                level = _determine_level([a21, a42, a63], pAdics, T, Integer(1) )
    else:
        if is_odd(n):
            k = ZZ((n+Integer(1) )/Integer(2) )
            a21 = S(E.a2()/pi)
            if char == Integer(2) :
                a6k = S(E.a6()/(pi**(Integer(2) *k+Integer(1) )))
                level = _determine_level([a21, a6k], pAdics, T, Integer(1) )
            else:
                a4k = S(E.a4()/(pi**(k+Integer(1) )))
                level = _determine_level([a21, a4k], pAdics, T, Integer(1) )
        else:
            k = ZZ((n+Integer(2) )/Integer(2) )
            if char == Integer(2) :
                a6k = S(E.a6()/(pi**(Integer(2) *k)))
                level = _determine_level([a6k], pAdics, T, Integer(1) )
            else:
                a3k = S(E.a3()/(pi**k))
                level = _determine_level([a3k], pAdics, T, Integer(1) )

    # Determine the transformations necessary
    if verbose > Integer(0) :
        print("Determining necessary transformations for " +
               str(T.count_children_at_level(level)) + " cases")
    for node in T.children_at_level(level):
        if n == Integer(1) :
            a42_ev = F(a42(node.representative()))
            if char == Integer(2) :
                change = a42_ev**sqrtPower
            else:
                a21_ev = F(a21(node.representative()))
                if char == Integer(3) :
                    change = -a42_ev / (Integer(2)  * a21_ev)
                else:
                    a63_ev = F(a63(node.representative()))
                    change = (a21_ev*a42_ev - Integer(9) *a63_ev)/(Integer(2) *(Integer(3) *a42_ev - a21_ev**Integer(2) ))
        else:
            if is_odd(n):
                a21_ev = F(a21(node.representative()))
                if char == Integer(2) :
                    a6k_ev = F(a6k(node.representative()))
                    square = a6k_ev / a21_ev
                    change = square**sqrtPower
                else:
                    a4k_ev = F(a4k(node.representative()))
                    change = - a4k_ev / (Integer(2)  * a21_ev)
            else:
                if char == Integer(2) :
                    a6k_ev = F(a6k(node.representative()))
                    change = a6k_ev**sqrtPower
                else:
                    a3k_ev = F(a3k(node.representative()))
                    change = - a3k_ev / F(Integer(2) )
                    
        if change in changeDict:
            changeDict[change].merge(node, from_root=True)
        else:
            Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
            Tn.merge(node, from_root=True)
            changeDict[change] = Tn

    # Performing the necessary transformations
    for (change, Tn) in changeDict.items():
        if n==Integer(1) :
            r = pi * F.lift(change)
            t = Integer(0) 
        elif is_odd(n):
            r = pi**k * F.lift(change)
            t = Integer(0) 
        else:
            r = Integer(0) 
            t = pi**k * F.lift(change)
        En = E.rst_transform(r, Integer(0) , t)
        newCase = case.copy()
        newCase.update(dict(next=f"Step7({n})", T=Tn, E=En,
                            urst=_urst_combine(urst, (Integer(1) , r, Integer(0) , t))))
        if verbose > Integer(0) :
            newCase['number'] = _get_case_number()
            print("Added case " + str(newCase['number']) + ":")
            print("Transform curve into")
            print(En)
            print("for")
            nT = Tn.minimum_full_level()
            resultstr = ""
            for N in Tn.children_at_level(nT):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(Tn.pAdics().prime()) +
                          "^" + str(nT))
            print(resultstr)
        result.put(newCase)
        
def _tate_step7sub(n, case, result, restrictions, **kwds):
    r"""Perform a step in the subalgorithm of step 7 of Tate's algorithm.

    In step $n$ of this subalgorithm we have different stopping
    conditions depending on the parity of $n$.

    If $n$ is odd we stop if the polynomial .. MATH::

        X^2 + \frac{a_3}{\pi^{(n+3)/2}}X - \frac{a_6}{\pi^{n+3}}

    has distinct roots over the algebraic closure of the residue field
    of the p-adics, i.e. if the valuation of $a_3^2 + 4 a_6$ is less
    than $n + 4$.

    If $n$ is even we stop if the polynomial .. MATH::

        \frac{a_2}{\pi} X^2 + \frac{a_4}{\pi^{(n+4)/2} X
        + \frac{a_6}{\pi^(n+3)}
    
    has distinct roots over the algebraic closure of the residue field
    of the p-adics, i.e. if the valuation of $a_4^2 - 4*a_2*a_6$ is
    less than $n + 5$.

    In both cases $a_2$, $a_3$, $a_4$ and $a_6$ are the appropiate
    invariants of the elliptic curve and $\pi$ is the uniformizer of
    the p-adics. If we don't stop in a step we continue to the next
    step of the subalgorithm.

    In case we stop in step $n$ of the subalgorithm we have:
    
    - Kodaira symbol: $I_n$
    
    - number of components: n + 5

    INPUT:

    - ``n`` -- The number of the step of the subalgorithm that will be
      performed

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    if (n > Integer(4) *pAdics.valuation(Integer(2) ) and
        len(restrictions) > Integer(0) ):
        newCase = case.copy()
        del newCase['next']
        newCase.update(dict(KS="In*"))
        _get_cases_invariant(S(E.c4()), pAdics, T, 'f', newCase,
                             result, **kwds)
    else:
        case_big = case.copy()
        case_big['next'] = f"Step7({n+1})T"
        case_small = case.copy()
        del case_small['next']
        case_small.update(dict(KS=f"I{n}*", m=Integer(5) +n))
        if is_odd(n):
            _get_two_cases_invariant(S(E.a3()**Integer(2)  + Integer(4) *E.a6()), pAdics,
                                     T, n + Integer(4) , case_big,
                                     case_small, result,
                                     **kwds)
        else:
            _get_two_cases_invariant(S(E.a4()**Integer(2)  - Integer(4) *E.a2()*E.a6()),
                                     pAdics, T, n + Integer(5) ,
                                     case_big, case_small,
                                     result, **kwds)
     
def _tate_step8_t(case, result, verbose=False):
    r"""Perform the transformation necessary before step 8 of Tate's
    algorithm.

    Transform the curve such that it has a Weierstrass equation of the
    form .. MATH::

        Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

    for which $a_2$ has valuation at least 2, $a_4$ has valuation at
    least 3 and $a_6$ has valuation at least 4.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step
    
    - ``verbose`` -- Verbosity argument

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    urst = case['urst']
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    pi = pAdics.uniformizer()
    
    if F.characteristic() != Integer(3) :
        # If the characteristic is not 3,
        # the transformation is always the same
        r = -E.a2() * F.lift(F(Integer(3) )**(-Integer(1) ))
        En = E.rst_transform(r, Integer(0) , Integer(0) )
        newCase = case.copy()
        newCase.update(dict(next='Step8', E=En,
                            urst=_urst_combine(urst, (Integer(1) , r, Integer(0) , Integer(0) ))))
        if verbose > Integer(0) :
            newCase['number'] = _get_case_number()
            print("Added case " + str(newCase['number']) + ":")
            print("Transform curve into")
            print(En)
            print("for")
            nT = T.minimum_full_level()
            resultstr = ""
            for N in T.children_at_level(nT):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(T.pAdics().prime()) +
                          "^" + str(nT))
            print(resultstr)
        result.put(newCase)
    else:
        # In characteristic 3
        # Find the integer s such that a^s is the cube root of a
        if R == ZZ and P in ZZ.ideal_monoid():
            cubertPower = ((ZZ.quotient(P.gen() - Integer(1) )(Integer(3) ))**(-Integer(1) )).lift()
        else:
            cubertPower = ((ZZ.quotient(ZZ(norm(P) - Integer(1) ))(Integer(3) ))**(-Integer(1) )).lift()

        # Determine the necessary transformations
        changeDict = dict()
        a63 = S(E.a6()/(pi**Integer(3) ))
        level = _determine_level([a63], pAdics, T, Integer(1) )
        if verbose > Integer(0) :
            print("Determining necessary transformation for " +
                   str(T.count_children_at_level(level)) + " cases")
        for node in T.children_at_level(level):
            cube = F(a63(node.representative()))
            change = cube**cubertPower
            if change in changeDict:
                changeDict[change].merge(node, from_root=True)
            else:
                Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
                Tn.merge(node, from_root=True)
                changeDict[change] = Tn

        # Performing the necessary transformations
        for (change, Tn) in changeDict.items():
            r = -pi * F.lift(change)
            En = E.rst_transform(r, Integer(0) , Integer(0) )
            newCase = case.copy()
            newCase.update(dict(next='Step8', T=Tn, E=En,
                                urst=_urst_combine(urst, (Integer(1) , r, Integer(0) , Integer(0) ))))
            if verbose > Integer(0) :
                newCase['number'] = _get_case_number()
                print("Added case " + str(newCase['number']) + ":")
                print("Transform curve into")
                print(En)
                print("for")
                nT = Tn.minimum_full_level()
                resultstr = ""
                for N in Tn.children_at_level(nT):
                    resultstr += str(N.representative()) + " "
                resultstr += ("modulo " + str(Tn.pAdics().prime()) +
                              "^" + str(nT))
                print(resultstr)
            result.put(newCase)
            
def _tate_step8(case, result, **kwds):
    r"""Perform step 8 of Tate's algorithm.

    Stop if the polynomial .. MATH::
    
        X^2 + \frac{a_3}{\pi^2} X + \frac{a_6}{\pi^4}

    has distinct roots in the algebraic closure of the residue field
    of the p-adics. Equivalently we stop if the valuation of $a_3^2 -
    4 a_6$ is less than 5. Here $a_3$ and $a_6$ are the
    respective invariants of the elliptic curve and $\pi$ is the
    uniformizer of the p-adics. Continue to step 9 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $IV^*$
    
    - number of components: 7

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    case_big = case.copy()
    case_big['next'] = 'Step9T'
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(KS="IV*", m=Integer(7) ))
    char = pAdics.characteristic()
    if char != Integer(3) :
        case_small['vDelta'] = Integer(8) 
        case_small['f'] = Integer(2) 
    _get_two_cases_invariant(S( ( E.a3() )**Integer(2)  + Integer(4)  * E.a6() ), pAdics,
                             T, Integer(5) , case_big, case_small, result,
                             **kwds)
        
def _tate_step9_t(case, result, verbose=False):
    r"""Perform the transformation necessary before step 9 of Tate's
    algorithm.

    Transform the curve such that it has a Weierstrass equation of the
    form .. MATH::

        Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

    for which $a_3$ has valuation at least 3 and $a_6$ has valuation
    at least 5.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step
    
    - ``verbose`` -- Verbosity argument

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    urst = case['urst']
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    pi = pAdics.uniformizer()
    
    if F.characteristic() != Integer(2) :
        # If the characteristic is not 2
        # the transformation is always the same
        t = -E.a3() * F.lift(F(Integer(2) )**(-Integer(1) ))
        En = E.rst_transform(Integer(0) , Integer(0) , t)
        newCase = case.copy()
        newCase.update(dict(next='Step9', E=En,
                            urst=_urst_combine(urst, (Integer(1) , Integer(0) , Integer(0) , t))))
        if verbose > Integer(0) :
            newCase['number'] = _get_case_number()
            print("Added case " + str(newCase['number']) + ":")
            print("Transform curve into")
            print(En)
            print("for")
            nT = T.minimum_full_level()
            resultstr = ""
            for N in T.children_at_level(nT):
                resultstr += str(N.representative()) + " "
            resultstr += ("modulo " + str(T.pAdics().prime()) +
                          "^" + str(nT))
            print(resultstr)
        result.put(newCase)
    else:
        # In characteristic 2
        # Determine the integers s such that a^s is the square root of a
        if R == ZZ and P in ZZ.ideal_monoid():
            sqrtPower = ( ( ZZ.quotient( P.gen() - Integer(1)  )(Integer(2) ) )**(-Integer(1) ) ).lift()
        else:
            sqrtPower = ( ( ZZ.quotient( ZZ(norm(P) - Integer(1) ) )(Integer(2) ) )**(-Integer(1) ) ).lift()
        if sqrtPower == Integer(0) :
            sqrtPower = Integer(1) 

        # Determine all necessary transformations
        changeDict = dict()
        a64 = S(E.a6()/(pi**Integer(4) ))
        level = _determine_level([a64], pAdics, T, Integer(1) )
        if verbose > Integer(0) :
            print("Determining necessary transformation for " +
                   str(T.count_children_at_level(level)) + " cases")
        for node in T.children_at_level(level):
            square = -F(a64(node.representative()))
            change = square**sqrtPower
            if change in changeDict:
                changeDict[change].merge(node, from_root=True)
            else:
                Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
                Tn.merge(node, from_root=True)
                changeDict[change] = Tn

        # Performing all necessary transformations
        for (change, Tn) in changeDict.items():
            t = -pi**Integer(2)  * F.lift(change)
            En = E.rst_transform(Integer(0) , Integer(0) , t)
            newCase = case.copy()
            newCase.update(dict(next='Step9', T=Tn, E=En,
                                urst=_urst_combine(urst, (Integer(1) , Integer(0) , Integer(0) , t))))
            if verbose > Integer(0) :
                newCase['number'] = _get_case_number()
                print("Added case " + str(newCase['number']) + ":")
                print("Transform curve into")
                print(En)
                print("for")
                nT = Tn.minimum_full_level()
                resultstr = ""
                for N in Tn.children_at_level(nT):
                    resultstr += str(N.representative()) + " "
                resultstr += ("modulo " + str(Tn.pAdics().prime()) +
                              "^" + str(nT))
                print(resultstr)
            result.put(newCase)
            
def _tate_step9(case, result, **kwds):
    r"""Perform step 9 of Tate's algorithm.

    Stop if the valuation of the invariant $a_4$ of the elliptic curve
    is less than 4. Continue to step 10 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $III^*$
    
    - number of components: 8

    - order of the group of components: 2

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    case_big = case.copy()
    case_big['next'] = 'Step10'
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(KS="III*", m=Integer(8) , c=Integer(2) ))
    char = pAdics.characteristic()
    if char != Integer(2) :
        case_small['vDelta'] = Integer(9) 
        case_small['f'] = Integer(2) 
    _get_two_cases_invariant(S(E.a4()), pAdics, T, Integer(4) , case_big,
                             case_small, result, **kwds)
            
def _tate_step10(case, result, **kwds):
    r"""Perform step 10 of Tate's algorithm.

    Stop if the valuation of the invariant $a_6$ of the elliptic curve
    is less than 6. Continue to step 11 otherwise.

    In case we stop we have:
    
    - Kodaira symbol: $II^*$
    
    - number of components: 9

    - order of the group of components: 1

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    case_big = case.copy()
    case_big['next'] = 'Step11'
    case_small = case.copy()
    del case_small['next']
    case_small.update(dict(KS="II*", m=Integer(9) , c=Integer(1) ))
    char = pAdics.characteristic()
    if char != Integer(2)  and char != Integer(3) :
        case_small['vDelta'] = Integer(10) 
        case_small['f'] = Integer(2) 
    _get_two_cases_invariant(S(E.a6()), pAdics, T, Integer(6) , case_big,
                             case_small, result, **kwds)
                               
def _tate_step11(case, result, verbose=False):
    r"""Perform step 11 of Tate's algorithm.

    Rescale the curve and start over at step 1.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step
    
    - ``verbose`` -- Verbosity argument

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    urst = case['urst']
    pi = pAdics.uniformizer()
    a1 = S(E.a1()/pi)
    a2 = S(E.a2()/(pi**Integer(2) ))
    a3 = S(E.a3()/(pi**Integer(3) ))
    a4 = S(E.a4()/(pi**Integer(4) ))
    a6 = S(E.a6()/(pi**Integer(6) ))
    E = EllipticCurve([a1,a2,a3,a4,a6])
    urst=_urst_combine(urst, (pi, Integer(0) , Integer(0) , Integer(0) ))
    newCase = case.copy()
    newCase.update(dict(next='Step1', E=E, E0=E, urst=urst,
                        urst0=urst))
    if verbose > Integer(0) :
        newCase['number'] = _get_case_number()
        print("Added case " + str(newCase['number']) + ":")
        print("Transform curve into")
        print(E)
        print("for")
        nT = T.minimum_full_level()
        resultstr = ""
        for N in T.children_at_level(nT):
            resultstr += str(N.representative()) + " "
        resultstr += ("modulo " + str(T.pAdics().prime()) +
                      "^" + str(nT))
        print(resultstr)
    result.put(newCase)
    
def _should_calculate_vDelta(case, restrictions):
    r"""Determine whether one should calculate the valuation of the
    discriminant.

    INPUT:

    - ``case`` -- The case for which this should be determined

    - ``restrictions`` -- Any restrictions on what should be computed

    OUTPUT:

    True if it is necessary to compute the valuation of the
    discriminant. False otherwise.

    """
    return (not ('vDelta' in case) and
            (len(restrictions) == Integer(0)  or 'discriminant' in restrictions or
             _should_calculate_m(case, restrictions) or
             _should_calculate_f(case, restrictions) or
             _should_calculate_n(case, restrictions)))

def _tate_calculate_vDelta(case, result, **kwds):
    r"""Compute the valuation of the discriminant.

    Compute the valuation of the discriminant of the elliptic curve.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    E = case['E']
    T = case['T']
    S = case['S']
    pAdics = case['pAdics']
    _get_cases_invariant(S(E.discriminant()), pAdics, T,
                         'vDelta', case, result, **kwds)

def _should_calculate_m(case, restrictions):
    r"""Determine whether one should calculate the number of components.

    INPUT:

    - ``case`` -- The case for which this should be determined

    - ``restrictions`` -- Any restrictions on what should be computed

    OUTPUT:

    True if it is necessary to compute the number of components. False
    otherwise.

    """
    return (not ('m' in case) and
            (len(restrictions) == Integer(0)  or
             _should_calculate_f(case, restrictions)))
    
def _tate_calculate_m(case, result):
    r"""Compute the number of components.

    Compute the number of components of the elliptic curve over the
    algebraic closure of the residue field, counted without
    multiplicity.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    """
    KS = case['KS']
    if KS.endswith('*'):
        if KS.startswith('IV'): #IV*
            case['m'] = Integer(7) 
        elif KS.startswith('III'): #III*
            case['m'] = Integer(8) 
        elif KS.startswith('II'): #II*
            case['m'] = Integer(9) 
        elif KS.startswith('I0'): #I0*
            case['m'] = Integer(5) 
        else: #In*
            case['m'] = case['vDelta'] - case['f'] + Integer(1) 
    else:
        if KS.startswith('IV'): #IV
            case['m'] = Integer(3) 
        elif KS.startswith('III'): #III
            case['m'] = Integer(2) 
        elif KS.startswith('II'): #II
            case['m'] = Integer(1) 
        elif KS.startswith('I0'): #II
            case['m'] = Integer(1) 
        else: #In
            case['m'] = case['vDelta']
    result.put(case)

def _should_calculate_c(case, restrictions):
    r"""Determine whether one should calculate the order of the group of
    components.

    INPUT:

    - ``case`` -- The case for which this should be determined

    - ``restrictions`` -- Any restrictions on what should be computed

    OUTPUT:

    True if it is necessary to compute the order of the group of
    components. False otherwise.

    """
    return (not ('c' in case) and len(restrictions) == Integer(0) )
    
def _tate_calculate_c(case, result, **kwds):
    r"""Compute the order of the group of components.

    Compute the number of components of the special fiber of the
    elliptic curve that are defined over the residue field of the
    p-adics.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`
    
    - ``verbose`` -- Verbosity argument

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    E = case['E']
    T = case['T']
    S = case['S']
    KS = case['KS']
    pAdics = case['pAdics']
    pi = pAdics.uniformizer()
    if KS.endswith('*'):
        if KS.startswith('IV'): #IV*
            if 'roots' in case:
                case['c'] = case['roots'] + Integer(1) 
                result.put(case)
            else:
                _get_number_of_roots_cases([S(Integer(1) ), S(E.a3()/(pi**Integer(2) )),
                                            S(-E.a6()/(pi**Integer(4) ))],
                                           pAdics, T, 'roots', case,
                                           result, **kwds)
        elif KS.startswith('III'): #III*
            case['c'] = Integer(2) 
            result.put(case)
        elif KS.startswith('II'): #II*
            case['c'] = Integer(1) 
            result.put(case)
        elif KS.startswith('I0'): #I0*
            if 'roots' in case:
                case['c'] = case['roots'] + Integer(1) 
                result.put(case)
            else:
                _get_number_of_roots_cases([S(Integer(1) ), S(E.a2()/pi),
                                            S(E.a4()/(pi**Integer(2) )),
                                            S(E.a6()/(pi**Integer(3) ))], pAdics,
                                           T, 'roots', case, result,
                                           **kwds)
        else: #In*
            if 'roots' in case:
                case['c'] = case['roots'] + Integer(2) 
                result.put(case)
            else:
                n = Integer(KS[Integer(1) :-Integer(1) ])
                if is_odd(n):
                    k = ZZ((n+Integer(3) )/Integer(2) )
                    _get_number_of_roots_cases([S(Integer(1) ),
                                                S(E.a3()/(pi**(k))),
                                                S(-E.a6()/(pi**(Integer(2) *k)))],
                                               pAdics, T, 'roots',
                                               case, result, **kwds)
                else:
                    k = ZZ((n+Integer(2) )/Integer(2) )
                    _get_number_of_roots_cases([S(E.a2()/pi),
                                                S(E.a4()/(pi**(k+Integer(1) ))),
                                                S(E.a6()/(pi**(Integer(2) *k+Integer(1) )))],
                                               pAdics, T, 'roots',
                                               case, result, **kwds)
    else:
        if KS.startswith('IV'): #IV
            if 'roots' in case:
                case['c'] = case['roots'] + Integer(1) 
                result.put(case)
            else:
                _get_number_of_roots_cases([S(Integer(1) ), S(E.a3()/pi),
                                            S(-E.a6()/(pi**Integer(2) ))],
                                           pAdics, T, 'roots', case,
                                           result, **kwds)
        elif KS.startswith('III'): #III
            case['c'] = Integer(2) 
            result.put(case)
        elif KS.startswith('II'): #II
            case['c'] = Integer(1) 
            result.put(case)
        elif KS.startswith('I0'): #II
            case['c'] = Integer(1) 
            result.put(case)
        else: #In
            if case['split']:
                case['c'] = case['vDelta']
                result.put(case)
            else:
                case['c'] = (Integer(1)  if is_odd(case['vDelta']) else Integer(2) )
                result.put(case)

def _should_calculate_f(case, restrictions):
    r"""Determine whether one should calculate the conductor exponent.

    INPUT:

    - ``case`` -- The case for which this should be determined

    - ``restrictions`` -- Any restrictions on what should be computed

    OUTPUT:

    True if it is necessary to compute the conductor exponent. False
    otherwise.

    """
    return (not ('f' in case) and (len(restrictions) == Integer(0)  or
                                       'conductor' in restrictions or
                                       'reduction_type' in restrictions))
    
def _tate_calculate_f(case, result):
    r"""Compute the conductor exponent.

    Compute exponent of the conductor of the elliptic curve at the
    prime of the p-adics.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    """
    case['f'] = case['vDelta'] - case['m'] + Integer(1) 
    result.put(case)

def _should_calculate_n(case, restrictions):
    r"""Determine whether one should calculate the number $n$ in the
    Kodaira symbol.

    INPUT:

    - ``case`` -- The case for which this should be determined

    - ``restrictions`` -- Any restrictions on what should be computed

    OUTPUT:

    True if it is necessary to compute the number $n$ in the Kodaira
    symbol. False otherwise.

    """
    return (case['KS'].count('n') > Integer(0)  and (len(restrictions) == Integer(0)  or
                                           'type' in restrictions))
    
def _tate_calculate_n(case, result):
    r"""Compute the number $n$ in $I_n$ or $I_n^*$.

    Compute the number $n$ if the elliptic curve at the prime of the
    p-adics has Kodaira symbol $I_n$ or $I_n^*$.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    """
    KS = case['KS']
    if KS.endswith('*'): #In*
        case['KS'] = KS.replace('n', str(case['vDelta'] - Integer(4)  - case['f']))
    else: #In
        case['KS'] = KS.replace('n', str(case['vDelta']))
    result.put(case)

def _should_calculate_split(case, restrictions):
    r"""Determine whether one should calculate whether the reduction is
    split.

    INPUT:

    - ``case`` -- The case for which this should be determined

    - ``restrictions`` -- Any restrictions on what should be computed

    OUTPUT:

    True if it is necessary to compute whether the reduction is
    split. False otherwise.

    """
    return ((_should_calculate_c(case, restrictions) or
             'reduction_type' in restrictions) and
            (case['f'] == Integer(1)  and not ('split' in case)))
            
def _tate_calculate_split(case, result, **kwds):
    r"""Compute if split or non-split multiplicative reduction.

    INPUT:

    - ``case`` -- A dictionary containing the data of the case on
      which to perform this step

    - ``result`` -- A queue in which to put the different cases of
      this step

    - ``variables`` -- A list of the variables of the polynomial or
      None if it should be determined from `S`

    - ``result`` -- A list of resulting cases to which the cases
      determined by this function should be added.

    - ``precision_cap`` -- A cap on the precision used for the
      variables

    """
    if 'roots' in case:
        if case['roots'] > Integer(0) :
            case['split'] = True
        else:
            case['split'] = False
        result.put(case)
    else:
        E = case['E']
        T = case['T']
        S = case['S']
        pAdics = case['pAdics']
        _get_number_of_roots_cases([S(Integer(1) ), S(E.a1()), S(-E.a2())],
                                   pAdics, T, 'roots', case, result,
                                   **kwds)
    
def _tate_finish(case, restrictions, result=[], variables=None,
                 input_data=False):
    r"""Turn a case into a final result

    INPUT:

    - ``case`` -- The case to be finalized
    
    - ``restrictions`` -- The restrictions on what should be computed

    - ``result`` -- A list to which the final result should be appended

    - ``variables`` -- The variables in the coefficients of the
      elliptic curve.

    OUTPUT:

    The given list 'result' to which we appended a tuple consisting of
    the information needed in the final result and a pAdicTree object
    that contains the values for the variables necessary to get this
    result.

    """
    tree = pAdicTree(variables=variables, root=case['T'])
    if len(restrictions) == Integer(0) :
        f = case['f']
        if f == Integer(0) :
            red_type = None
        elif f == Integer(1) :
            if case['split']:
                red_type = Integer(1) 
            else:
                red_type = -Integer(1) 
        else:
            red_type = Integer(0) 
        R = case['E0'].base_ring()
        urst = tuple(R(x) for x in case['urst0'])
        urst_inv = tuple(R(x) for x in _urst_invert(case['urst0']))
        myresult = FreyCurveLocalData(case['E0'],
                                      case['T'].pAdics().prime_ideal(),
                                      f, case['vDelta'],
                                      KodairaSymbol(case['KS']),
                                      case['c'], red_type, urst,
                                      urst_inv)
    else:
        myresult = []
        for r in restrictions:
            if r == 'conductor':
                myresult.append(case['f'])
            if r == 'reduction_type':
                f = case['f']
                if f == Integer(0) :
                    myresult.append(None)
                elif f == Integer(1) :
                    if case['split']:
                        myresult.append(Integer(1) )
                    else:
                        myresult.append(-Integer(1) )
                else:
                    myresult.append(Integer(0) )
            if r == 'discriminant':
                myresult.append(case['vDelta'])
            if r == 'type':
                myresult.append(KodairaSymbol(case['KS']))
            if r == 'minimal_model':
                myresult.append(case['E0'])
            if r == 'isomorphism':
                R = case['E0'].base_ring()
                urst = tuple(R(x) for x in case['urst0'])
                urst_inv = tuple(R(x) for x in
                                 _urst_invert(case['urst0']))
                myresult.append((urst, urst_inv))
    if input_data:
        my_result = ((case['Estart'], case['pAdics']), my_result)
    result.append((myresult, tree))
    return result
           
def _tate_cleanup(cases):
    r"""Cleanup the final result to return

    INPUT:

    - ``cases`` -- A list of pairs consisting of possible return
      values and a pAdicTree containing the values of the variables
      necessary to get those return values.

    OUTPUT:

    The result of the function as a ConditionalValue

    """
    result = []
    for case in cases:
        flag = True
        for case0 in result:
            if case0[Integer(0) ] == case[Integer(0) ]:
                # We can merge here to prevent unneccesary copying.
                # Since we would discard both old trees anyways this
                # is not a problem.
                case0[Integer(1) ].pAdic_tree()._root.merge(case[Integer(1) ]._root)
                flag = False
                break
        if flag:
            result.append((case[Integer(0) ], TreeCondition(case[Integer(1) ])))
    return ConditionalValue(result)

def _check_elliptic_curve(E):
   r"""Check whether input is an elliptic curve"""
   if not isinstance(E, EllipticCurve_generic):
       raise ValueError("%s is not an elliptic curve."%(E,))

def _init_coefficient_ring(coefficient_ring, E):
    r"""Initialize the coefficient ring

    INPUT:

    - ``coefficient_ring`` -- The argument coefficient ring

    - ``E`` -- The elliptic curve

    OUTPUT:

    The coefficient ring to be used

    """
    if coefficient_ring is None:
        coefficient_ring = E.base_ring()
    for a in E.a_invariants():
        if a not in coefficient_ring:
            raise ValueError("%s is not part of %s."%(a, coefficient_ring))
    return coefficient_ring

def _init_pAdics(pAdics, ring, prime, coefficient_ring):
    if pAdics is None:
        if ring is None:
            ring = coefficient_ring.base_ring()
        if prime is None:
            raise ValueError("At least the argument pAdics or " +
                             "the argument prime should be given.")
        pAdics = pAdicBase(ring, prime)
    if not isinstance(pAdics, pAdicBase):
        raise ValueError("%s is not a pAdicBase object."%(pAdics,))
    return pAdics
 
def _init_polynomial_ring(coefficient_ring, pAdics):
    return coefficient_ring.change_ring(pAdics.number_field())
    
def _init_variables_tate(polynomial_ring):
    return list(polynomial_ring.gens())
    
def _init_initial_values(initial_values, pAdics, variables):
    if initial_values is None:
        initial_values = pAdicTree(pAdics=pAdics, full=True,
                                   variables=variables)
    if not isinstance(initial_values, pAdicTree):
        raise ValueError("%s is not a p-adic tree."%(initial_value,))
    if not pAdics.is_extension_of(initial_values.pAdics()):
        raise ValueError(str(pAdics) + " does not extend the p-adics of " +
                         str(initial_values))
    return initial_values
    
def _init_case(T, E, S, pAdics, verbose):
    firstCase = dict(next='Step1', T=T.root(), E=E, E0=E,
                     pAdics=pAdics, S=S,
                     urst=(Integer(1) , Integer(0) , Integer(0) , Integer(0) ),
                     urst0=(Integer(1) , Integer(0) , Integer(0) , Integer(0) ))
    if verbose > Integer(1) :
        _init_case_number()
        firstCase['number'] = "0"
    return firstCase
    
def _init_str_list(str_list):
    if isinstance(str_list, str) or not hasattr(str_list, '__iter__'):
        str_list = [str_list]
    for s in str_list:
        if not isinstance(s, str):
            raise ValueError("%s is not a string."%(s,))
    return str_list

def _urst_combine(urst1, urst2):
    r"""Combine two changes of a Weierstrass model into one

    INPUT:

    - ``urst1`` -- A tuple describing the first change in Weierstrass
      model
    
    - ``urst2`` -- A tuple describing the second change in Weierstrass
      model

    OUTPUT:

    A tuple describing a single change in Weierstrass model that has
    the same result as first doing the change described by `urst1` and
    then the change described by `urst2`.

    """
    u1, r1, s1, t1 = urst1
    u2, r2, s2, t2 = urst2
    return (u1*u2, u1**Integer(2) *r2 + r1, u1*s2 + s1, u1**Integer(2) *r2*s1 + u1**Integer(3) *t2 + t1)

def _urst_invert(urst):
    r"""Invert the change of Weierstrass model

    INPUT:

    - ``urst`` -- A tuple describing a change in Weierstrass model
      from an elliptic curve `E1` to an elliptic curve `E2`

    OUTPUT:

    A tuple describing a change in Weierstrass model from `E2` to `E1`
    that is the inverse of `urst`.

    """
    u, r, s, t = urst
    uinv = u**(-Integer(1) )
    return (uinv, -uinv**Integer(2)  * r, -uinv * s, uinv**Integer(3) *(r*s - t))

