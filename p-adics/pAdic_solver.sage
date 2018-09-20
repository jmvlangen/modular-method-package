r"""
Methods for finding solutions to polynomial equations over local fields

Given a polynomial with coefficients in some number field
and a finite prime P, can compute the roots of this polynomial in
the completion at P up to a certain precision and return the result
as a p-adic tree.

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Joey van Langen (2018-07-13): initial version

"""

# ****************************************************************************
#       Copyright (C) 2018 Joey van Langen <j.m.van.langen@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import Infinity

from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial

def find_pAdicRoots(polynomial, pAdics=None, ring=None, prime=None,
                    variables=None, value_tree=None,
                    precision=20, verbose=False, give_list=False,
                    precision_cap=20):
    r"""
    Finds the pAdic roots of a polynomial.
    
    INPUT:
    
    - ``polynomial`` -- A polynomial (in several variables) of which we
      want to find pAdic roots.
      
    - ``pAdics`` -- A pAdicBase object (default: pAdicBase(ring, prime))
      which determines the pAdics to be used. Note that the arguments
      ring and prime will be ignored when this is given, and that
      the argument prime should be given is this argument was not.
      Furthermore the number field and prime of this pAdicBase object
      should satisfy the conditions given under ``ring`` and ``prime``
      respectively.
      
    - ``ring`` -- A ring (default: polynomial.base_ring()) over which
      we want to work p-adicly. Note that this should not be a p-adic
      ring, but rather a number field or ring of integers containing
      the coefficients of the polynomial and the prime over which we
      want to work p-adicly.
      
    - ``prime`` -- A prime ideal (default: None) of the argument
      ``ring``. If ring is a number field we mean by this a prime
      ideal of its ring of integers. This will be the prime over
      which we will work p-adicly.
      
    - ``variables`` -- An iterable object (default: polynomial.variables())
      containing the variables for which we want to substitute p-Adic
      values to obtain roots.
    
    - ``value_tree`` -- A pAdicSharedTree (default: A full pAdicSharedTree
      with width len(variables)) consisting of tuples that contain the
      pAdic values that the variables may attain. Note that the order
      of the values in these tuples must correspond to the order of the
      variables in the argument variables.
      
    - ``precision`` -- A non-negative integer (default: 20) determining
      to what precision p-Adic roots should be found. A precision of n
      means that the function finds all p-adic values that form a root
      modulo P^n where P is the prime ideal over which we do p-adics.
      
    - ``verbose`` -- A boolean value or an integer (default: False).
      When set to True or any value larger then zero will print comments
      to stdout about the computations being done whilst busy. If set
      to False or 0 will not print such comments. If set to any negative
      value will also prevent the printing of any warnings.
      If this method calls any method that accepts an argument verbose
      will pass this argument to it. If such a method fulfills a minor
      task within this method and the argument verbose was larger than
      0, will instead pass 1 less than the given argument. This makes it
      so a higher value will print more details about the computation
      than a lower one.
      
    - ``give_list`` -- A boolean value (default: False). When set to
      True the result of this function will be a list of trees. A
      tuple will be in the i-th tree in this list, if it is a root of
      the polynomial up to (precision - i) and is not a root for any
      smaller i >= 0.
      
    - ``precision_cap'' -- An integer (default:20). The highest precision
      for which to actually calculate. If this is set too high, the
      calculation might take a long time. Note that this is simply a cap
      on the size of the p-adic approximation of the variables, not of
      an approximation of the root itself as is precision!
      
    OUTPUT:
    
    A tuple consisting of
    
    - A pAdicSharedTree containing all the tuples of p-adic values that
      are roots of the given polynomial up to the given precision.
      
    - A pAdicSharedTree containing all the tuples of p-adic values that
      are not in the first tree, but that were given as values.
      
    If, however, the option give_list was set to True the output will
    be a tuple consisting of
    
    - A list in which each entry is a pAdicSharedTree. The n-th entry of
      this list contains the p-adic tuples that are roots up to precision
      n + min_prec (see second return value), but not up to any higher
      precision. The only exception to this rule are the roots of the
      precision given by the argument precision as the program does not
      check for roots above that level. Note that such an entry might
      not necessarily exist.
    - An integer, named min_prec, which gives the lowest precision in which
      roots can be found. Note that this value may be infinity, in which
      case the list might be empty.
    """
    _check_polynomial(polynomial)
    variables = _init_variables(variables, polynomial)
    pAdics = _init_pAdics(pAdics, ring, prime, polynomial)
    value_tree = _init_value_tree(value_tree, pAdics, variables)
    polynomial, least_power = _modify_polynomial(polynomial, pAdics, variables)
    multiplicity = pAdics.extension_multiplicity(value_tree.pAdics())
    precision, var_precision = _init_precision(precision, precision_cap,
                                               least_power, multiplicity,
                                               verbose)
    polynomial_derivatives = _init_derivatives(polynomial, variables)

    if verbose > 0:
        print "Finding roots of %s modulo %s^%s."(polynomial,
                                                  pAdics.prime_ideal()._repr_short(),
                                                  precision)
       
    result = _find_pAdicRoots(polynomial, polynomial_derivatives, value_tree,
                              precision, var_precision, multiplicity, pAdics,
                              verbose=(verbose-1 if verbose > 0 else verbose),
                              give_list=give_list,
                              quit_on_empty=give_list)
    if give_list:
        return result, least_power
    else:
        return result

def _pAdic_fill_in_check(f, T, prec, pAdics, verbose=False):
    if verbose > 0:
        print "Checking for roots at level 1"
    Tno = pAdicNode(pAdics=T.pAdics(), width=T.width)
    for node in T.children_at_level(1):
        if pAdics.valuation(f(node.representative())) < prec:
            Tno.merge(node, from_root=True)
            node.remove()
    return T, Tno

def _pAdic_convert_function(pAdics1, pAdics2, step_size):
    phi = pAdics1.extension_vector_space(pAdics2)
    F = pAdics1.residue_field()
    G = pAdics2.residue_field()
    def result(a):
        b = pAdics1.power_series(a, precision=step_size)
        d = []
        for c in b:
            d.extend(list(phi(c)))
        return (G^len(d))(d)
    return result
        
def _pAdic_hensel_check(f, fder, T, level, step_size, pAdics, verbose=False):
    if verbose > 0:
        print "Checking for roots at level %d"%level
    Tyes = pAdicNode(pAdics=T.pAdics(), width=T.width)
    phi = _pAdic_convert_function(pAdics, T.pAdics(), step_size)
    F = T.pAdics().residue_field()
    pi = T.pAdics().uniformizer()
    if verbose > 1:
        print "Start loop of %d cycles"%T.count_children_at_level(level-1)
    for node in T.children_at_level(level - 1):
        fa = phi(f(node.representative())/(pi^(level-1)))
        fdera = [phi(fderi(node.representative())) for fderi in fder]
        M = matrix(fdera)
        try:
            c = M.solve_left(-fa)
            for c0 in M.kernel():
                cfs = tuple([F.lift(ci) for ci in (c+c0)])
                if node.children.contains(cfs):
                    child = node.children.get(cfs)
                    Tyes.merge(child, from_root=True)
                    child.remove()
        except ValueError:
            pass
    return Tyes, T
    
def _pAdic_level_check(f, fder, T, level, step_size, pAdics, verbose=False):
    if level <= 0:
        return (T,)
    if level == 1:
        return _pAdic_fill_in_check(f, T, step_size, pAdics, verbose=verbose)
    if level > 1:
        return _pAdic_hensel_check(f, fder, T, level, step_size, pAdics,
                                      verbose=verbose)
    raise ValueError("The level %s is not valid."%level)
        
def _find_pAdicRoots(f, fder, T, prec, var_prec, m, pAdics, verbose=False,
                     give_list=False, quit_on_empty=False):
    r"""
    The recursive implementation of :func:find_pAdicRoots
    
    Note that this does not check or default arguments!
    Only for internal use.
    """
    if give_list:
        Tno = []
    else:
        Tno = pAdicNode(pAdics=T.pAdics(), width=T.width)
    if prec > 0:
        for level in range(1,var_prec+1):
            if T.is_empty():
                if quit_on_empty or not give_list:
                    if give_list:
                        return Tno
                    else:
                        return T, Tno
                level_result = (T,T)
            else:
                step_size = min(m, prec-m*(level-1))
                level_result = _pAdic_level_check(f, fder, T,
                                                  level, step_size, pAdics,
                                                  verbose=verbose)
            T = level_result[0]
            if len(level_result) > 1:
                if give_list:
                    Tno.append(level_result[1])
                else:
                    Tno.merge(level_result[1])
    if give_list:
        Tno.append(T)
        return Tno
    else:
        return T, Tno
        
def _check_polynomial(polynomial):
    if not isinstance(polynomial, Polynomial) and \
       not isinstance(polynomial, MPolynomial):
         raise ValueError("%s is not a polynomial."%(polynomial,))
    
def _init_variables(variables, polynomial):
    if variables is None:
        variables = polynomial.variables()
    if not hasattr(variables, '__iter__'):
        raise ValueError("%s is not a list of variables"%(variables,))
    for v in polynomial.variables():
        if v not in variables:
            raise ValueError("%s is a variable of the polynomial, but not listed as a variable"%(v,))
    return variables
    
def _init_pAdics(pAdics, ring, prime, polynomial):
    if pAdics is None:
        if ring is None:
            ring = polynomial.base_ring()
        if prime is None:
            raise ValueError("At least the argument pAdics or the argument prime should be given.")
        pAdics = pAdicBase(ring, prime)
    if not isinstance(pAdics, pAdicBase):
        raise ValueError("%s is not a pAdicBase object."%(pAdics,))
    return pAdics
    
def _init_value_tree(value_tree, pAdics, variables):
    if value_tree is None:
        value_tree = pAdicNode(pAdics=pAdics, full=True, width=len(variables))
    if isinstance(value_tree, pAdicTree):
        value_tree = value_tree.root().copy()
    if not isinstance(value_tree, pAdicNode):
        raise ValueError("%s does not define a valid p-adic tree."%(value_tree,))
    if not pAdics.is_extension_of(value_tree.pAdics()):
        raise ValueError("%s does not match the p-adics of %s."%(pAdics,
                                                                 value_tree))
    if value_tree.width != len(variables):
        raise ValueError("%s does not have as many entries as there are variables"%(value_tree))
    return value_tree

def _modify_polynomial(polynomial, pAdics, variables):
    S = PolynomialRing(pAdics.number_field(), variables)
    polynomial = S(polynomial)
    least_power = min([Infinity] +
                      [pAdics.valuation(c) for c in polynomial.coefficients()])
    if least_power < Infinity:
        polynomial = pAdics.uniformizer()^(-least_power) * polynomial
    return polynomial, least_power
           
def _init_precision(precision, precision_cap, least_power, multiplicity, verbose):
    if not precision in ZZ:
        raise ValueError("The precision should be an integer, not %s."%(precision,))
    if not precision_cap in ZZ:
        raise ValueError("The precision cap should be an integer, not %s"%(precision_cap,))
    precision = precision - least_power
    if precision > multiplicity * precision_cap and verbose >= 0:
        precision = multiplicity * precision_cap
        print "Warning: Lowering precision on root to %d to accomodate for precision_cap on variables"%(precision+least_power)
    if precision < 0:
        precision = -1
    return precision, ceil(precision/multiplicity)
    
def _init_derivatives(polynomial, variables):
    result = []
    for v in variables:
        v = polynomial.parent()(v)
        result.append(polynomial.derivative(v))
    return result   
        
