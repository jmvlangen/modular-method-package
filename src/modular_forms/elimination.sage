r"""Functions to eliminate newforms corresponding to Frey curves

Given a diophantine problem for which we can associate to a solution a
Frey–Hellegouarch curve, we can prove the non-existence of that
solution by showing that the mod $l$ galois representation of that
curve is the same as the mod $l$ galois representation of a newform in
some finite list and procedurally eliminate these newforms. This file
provides implementations of the often used methods to eliminate these
newforms.

For the most part elimination is done by comparing the traces of a
newform with the possible traces of the corresponding
Frey–Hellegouarch curve. The method :func:`eliminate_by_trace` does
this for the frobenius map above a single prime, whilst the method
:func:`eliminate_by_traces` does the same for multiple primes at a
time. Note that both function can also work with multiple Frey curves
at a time to improve elimination.

Besides standard comparison of traces, we can also use the method of
Kraus implemented in :meth:`kraus_method` to eliminate newforms for a
specific prime $l$.

Another method of eliminating newforms is provided by the function
:func:`eliminate_cm_forms` which eliminates all the newforms with
complex multiplication if the corresponding Frey curve does not have
complex multiplication.

Lastly there is a method $eliminate_primes$ which allows one to
eliminate newforms for which we can determine that a certain mod $l$
galois representation correspondence can not exist by some other
means.

Every elimination method in this file will return a uniform output
that can again be used again as a input for another elimination
method. This output is a list of tuples, where each tuple contains a
newform corresponding to each Frey curve on which the method was
applied and as a last entry an integer divisible by all primes $l$ for
which the mod $l$ of each newform in the tuple might still correspond
to the mod $l$ representation of the corresponding Frey curve.

The function :func:`combine_newforms` allows one to combine the
different outputs of elimination methods into a single output. This is
usefull if one wants to perform elimination methods first on several
Frey curves individually and then use the results to perform some
elimination methods in which you use all curves simultaneously.

EXAMPLES:

TODO

AUTHORS:

- Joey van Langen (2019-03-05): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import itertools
from functools import reduce

from sage.arith.functions import lcm
from sage.arith.misc import gcd
from sage.rings.finite_rings.integer_mod import mod
from sage.rings.fast_arith import prime_range
from sage.rings.ideal import is_Ideal

from sage.misc.cachefunc import cached_function
from sage.misc.misc_c import prod as product

from modular_method.padics.pAdic_base import pAdicBase

from modular_method.number_fields.field_constructors import common_embedding_field

from modular_method.elliptic_curves.Qcurves import Qcurve_base as Qcurve
from modular_method.elliptic_curves.frey_curves import FreyCurve

from modular_method.diophantine_equations.conditions import CongruenceCondition
from modular_method.diophantine_equations.conditions import TreeCondition
from modular_method.diophantine_equations.conditions import ConditionalValue
from modular_method.diophantine_equations.conditions import conditional_product
from modular_method.diophantine_equations.conditions import apply_to_conditional_value

from sage.all import ZZ, QQ, Integer, Integers

def _init_elimination_data(curves, newforms, condition):
    r"""Initialize the data used by the different elimination methods.

    INPUT:

    - ``curves`` -- The argument `curves`

    - ``newforms`` -- The argument `newforms`

    - ``condition`` -- The argument `condition`

    OUTPUT:

    - A tuple of Frey curves on which the elimination method should be
      applied.

    - A list of tuples, wherein each tuple consists of one newform for
      each Frey curve in the previous tuple and as last some integer
      divisible by all primes $l$ for which the mod $l$ galois
      representations of these newforms and their corresponding Frey
      curves might still agree.

    - A condition on the variables of the Frey curves given first.

    """
    if isinstance(curves, FreyCurve):
        curves = (curves,)
    else:
        curves = tuple(curve for curve in curves)
        parameters = None
        if len(curves) == 0:
            raise ValueError("At least one curve should be provided.")            
        for curve in curves:
            if not isinstance(curve, FreyCurve):
                raise ValueError("%s is not a Frey curve"%(curve,))
            if parameters is None:
                parameters = curve.parameters()
            elif parameters != curve.parameters():
                raise ValueError(str(curves[0]) + " and " + str(curve) +
                                 " don't have the same parameters")
    newforms = _init_newform_list(newforms, curves)
    if condition is None:
        condition = curves[0]._condition
    return curves, newforms, condition

@cached_function
def _init_traces(curves, condition, primes, powers, precision_cap, verbose):
    r"""Initialize the traces of Frobenius of some Frey curves.

    INPUT:

    - ``curves`` -- A tuple containing the Frey curves

    - ``condition`` -- A condition on the variables of the Frey curves
      that should hold.

    - ``primes`` -- A tuple of finite primes, one for each Frey curve,
      at which the traces should be determined.

    - ``powers`` -- A tuple of positive integers, one for each Frey
      curve, the powers of the frobenius elements to be considered.

    - ``precision_cap`` -- The maximal precision to be used on the variables.

    - ``verbose`` -- Verbosity argument

    OUTPUT:

    A list of tuples wherein each tuple contains a possible
    combination of the traces of Frobenius of the various Frey
    curves. Each i-th entry in such a tuple is a possible trace of
    Frobenius at the i-th prime for the i-th Frey curve.

    NOTE:

    This function is cached as higher level functions might request
    the same value of this function multiple times. If they do they
    should clear the cache after it is certain no similar calls to
    this function will be made.

    """
    traces = [curves[i].trace_of_frobenius(primes[i], power=powers[i],
                                           condition=condition,
                                           precision_cap=precision_cap,
                                           verbose=(verbose - 1 if verbose > 0
                                                    else verbose))
              for i in range(len(curves))]
    traces = conditional_product(*traces)
    traces = ConditionalValue([(val, con) for val, con in traces
                               if not con.never()])
                                       
def _init_newform_list(newforms, curves):
    """Initialize a list of newforms associated to given Frey curves

    INPUT:

    - ``newforms`` -- The argument `newforms` of an elimination
      function

    - ``curves`` -- A tuple or list of the Frey curves

    OUTPUT:

    A list of tuples, such that each tuple contains one newform
    corresponding to each Frey curve and as a last entry an
    integer. The integer is divisible by all those primes $l$ such
    that the mod $l$ galois representations of the newforms and those
    of the Frey curves might still agree.

    """
    if newforms is None and len(curves) == 1:
        newforms = apply_to_conditional_value(lambda x: list(x),
                                              curves[0].newform_candidates())
    if isinstance(newforms, tuple):
        if len(newforms) != len(curves):
            raise ValueError("Expected " + str(len(curves)) +
                             " lists of newforms, but got " +
                             str(len(newforms)))
        newforms = [(apply_to_conditional_value(lambda x: list(x),
                                                curves[i].newform_candidates())
                     if newforms[i] is None
                     else newforms[i])
                    for i in range(len(curves))]
        newforms = conditional_product(*newforms)
        newforms = apply_to_conditional_value(lambda nfs:
                                              itertools.product(*nfs),
                                              newforms)
    if isinstance(newforms, ConditionalValue):
        return apply_to_conditional_value(lambda nfs:
                                          _init_newform_list(nfs, curves),
                                          newforms)
    newforms = [nfs for nfs in newforms]
    for i in range(len(newforms)):
        if not isinstance(newforms[i], tuple):
            if len(curves) != 1:
                raise ValueError("Expected " + str(len(curves)) +
                                 " newforms per entry, but got 1")
            newforms[i] = (newforms[i], ZZ(0))
        if len(newforms[i]) == len(curves) + 1:
            pass
        elif len(newforms[i]) == len(curves):
            newforms[i] = tuple(list(newforms[i]) + [0])
        else:
            raise ValueError("Expected " + str(len(curves)) +
                             " newforms per entry, but got " +
                             str(len(newforms[i])))
    return newforms

def eliminate_by_trace(curves, newforms, prime, B=0, condition=None,
                       precision_cap=1, verbose=False):
    r"""Eliminate newforms associated to Frey curves by the trace of
    frobenius at a given prime.

    Let $E$ be a Frey curve, $f$ be a newform of which the mod $l$
    galois representations might agree with the mod $l$ galois
    representation of $E$, and $F$ be a frobenius element of the
    common domain of the afore mentioned representations. By comparing
    the traces of these representations at $F$ we can check if the
    representations are indeed the same. If not we eliminate the
    newform for mod $l$ representations. This method will perform this
    elimination for a single Frobenius element.

    Note that the value of the trace of Frobenius for the Frey curve
    might depend on the value of the variables. In this case we
    eliminate the newform if the trace of Frobenius for the newform
    can not correspond to any of the possible traces for the elliptic
    curve. If we have multiple Frey curves we compare traces
    simultaneously for each combination of associated newforms, in
    which case we eliminate a tuple of newforms if there is no
    associated possible tuple of traces.

    INPUT:

    - ``curves`` -- A Frey curve or a list/tuple of Frey curves that
      share the same parameters.

    - ``newforms`` -- A list of newforms that are associated to the
      given Frey curve, i.e. such that their mod-l representations
      could be isomorphic to the mod-l representation of the Frey
      curve.  Instead of a list one can also give the value None, in
      which case the list will be retrieved from the curve by the
      method :meth:`newform_candidates` of the given Frey curve.

      If multiple Frey curves are given this argument should be a
      tuple of length equal to the number of given Frey curves,
      containing for each Frey curve an entry as described above in
      the same order as the Frey curves are given.

      Alternatively the input may also be a list of tuples wherein
      each tuple has length equal to the number of Frey curves and
      each entry thereof is a newform associated to the corresponding
      Frey curve.

      In the last input format, each tuple may also be 1 longer than
      the number of Frey curves, in which case the last entry must be
      an integer divisible by all primes l for which the mod-l
      representations of the given elliptic curves may still be
      isomorphic to the corresponding newforms in this tuple. This can
      be used to chain multiple elimination methods.

      Everywhere where a list is expected, one may also give a
      ConditionalValue of which the possible values are lists of the
      expected format. The method will be applied to each entry
      individually and the condition will be replaced by the condition
      that both the given condition and the condition at which that
      value is attained hold.

    - ``prime`` -- A finite prime of the field in which the parameters
      of the Frey curves live. This will be the prime underlying the
      Frobenius element for which the traces of the galois
      representations at that element should be compared. It should be
      give as a prime number if the field of the parameters is the
      rationals and a prime ideal otherwise. Note that the Frobenius
      elements might be chosen for primes lying above these primes if
      the base fields of the Galois representations considered are
      larger than the field of the parameters.

    - ``B`` -- A non-negative integer (default: 0). Will only
      eliminate newforms based on their mod $l$ representation if $l$
      divides this number.

    - ``condition`` -- A Condition giving the restrictions on the
      parameters of the given Frey curve(s) that should be
      considered. By default this will be set to the condition stored
      in the first given Frey curve.

    - ``precision_cap`` -- A non-negative integer (default: 1) giving
      the maximal precision used to compute the reduction type of a
      Frey curve at a prime. Note that setting this value to a value
      higher than 1 might result in a long computation time and heavy
      memory usage and is strongly discouraged.

    - ``verbose`` -- A boolean value or an integer (default:
      False). When set to True or any value larger then zero will
      print comments to stdout about the computations being done
      whilst busy. If set to False or 0 will not print such comments.
      If set to any negative value will also prevent the printing of
      any warnings.  A higher value will cause more messages to be
      printed.

    OUTPUT:
    
    A list of tuples, wherein each tuple contains:

    - for each frey curve given a corresponding newform for which an
      isomorphism between the mod $l$ galois representations of that
      curve and that newform could not be excluded by this method.

    - as its last entry an integer that is divisible by all prime
      numbers $l$ for which the mentioned mod $l$ galois
      representations could still exist. Note that this number is
      always 0 if there is infinitely many such $l$ remaining.

    .. SEEALSO::

        :func:`eliminate_by_traces`,
        :func:`kraus_method`,
        :func:`eliminate_cm_forms`,
        :func:`eliminate_primes`

    """
    curves, newforms, condition = _init_elimination_data(curves, newforms,
                                                         condition)
    K = curves[0]._R.fraction_field()
    if not ((K == QQ and prime in ZZ and prime > 0 and prime.is_prime()) or
            (is_Ideal(prime) and prime in K.ideal_monoid() and prime.is_prime())):
        raise ValueError("Argument %s is not a valid prime."%(prime,))
    if not (B in ZZ and B >= 0):
        raise ValueError("Argument %s is not a non-negative integer."%(B,))
    if B != 0:
        B = B.prime_factors()
    return _eliminate_by_trace(curves, newforms, prime, B, condition,
                               precision_cap, verbose)

def _single_elimination(E, KE, LE, nfs, p, prime, pE, B, Bprod, C,
                        prec_cap, verbose):
    r"""Perform the elimination for a single tuple of newforms

    INPUT:

    - ``E`` -- The Frey curves

    - ``KE`` -- The base fields of the Galois representations
      associated to the Frey curves

    - ``LE`` -- The coefficient fields of the Galois representations
      associated to the Frey curves

    - ``nfs`` -- A tuple of newforms, one for each Frey curve, and an
      integer

    - ``p`` -- The prime number below `prime`

    - ``prime`` -- The base prime over which to do elimination

    - ``pE`` -- A tuple of primes above `prime`, one for each Frey
      curve in the respective base field of the Galois representation

    - ``B`` -- Zero or a list of prime numbers, the primes to be
      eliminated.

    - ``Bprod`` -- The product of the elements of `B`

    - ``C`` -- The condition to be used

    - ``prec_cap`` -- The precision cap on the variables

    - ``verbose`` -- Verbosity argument

    """
    nE = len(E)
    Bold = nfs[-1]
    if (gcd(Bprod, Bold) == 1 or (B == 0 or Bold != 0)):
        # Nothing to do
        return Bold
    apf = [nfs[i].trace_of_frobenius(p, power=powers[i])
           for i in range(nE)]
    Kf = [nfs[i].base_field() for i in range(nE)]
    Kcom = [common_embedding_field(KE[i], Kf[i], give_maps=True)
            for i in range(nE)]
    KphiE = [Kcom[i][1] for i in range(nE)]
    Kphif = [Kcom[i][2] for i in range(nE)]
    Kcom = [Kcom[i][0] for i in range(nE)]
    pcom = [(QQ(pE[i]) if Kcom[i] == QQ
             else Kcom[i].prime_above(KphiE[i](pE[i])))
            for i in range(nE)]
    pf = [(p if Kf[i] == QQ
           else next(P for P in Kf[i].primes_above(pf[i])
                     if Kcom[i].prime_above(Kphif[i](P)) == pcom[i]))
           for i in range(nE)]
    if any(pf[i].divides(nfs[i].level()) for i in range(nE)):
        # Can not do this prime, so skip
        return Bold
    ramdegE = [(1 if Kcom[i] == QQ else KphiE[i](pE[i]).valuation(pcom[i]))
               for i in range(nE)]
    ramdegf = [(1 if Kcom[i] == QQ else Kphif[i](pf[i]).valuation(pcom[i]))
               for i in range(nE)]
    resdegE = [(1 if Kcom[i] == QQ else
                (pcom[i].residue_class_degree() if KE[i] == QQ else
                 pcom[i].residue_class_degree() / pE[i].residue_class_degree()))
               for i in range(nE)]
    resdegf = [(1 if Kcom[i] == QQ else
                (pcom[i].residue_class_degree() if Kf[i] == QQ else
                 pcom[i].residue_class_degree() / pf[i].residue_class_degree()))
               for i in range(nE)]
    powcom = [lcm(ramdegE[i], ramdegf[i]) for i in range(nE)]
    powE = [resdegE[i]*powcom[i] for i in range(nE)]
    powf = [resdegE[i]*powcom[i] for i in range(nE)]
    apE_ls = _init_traces(E, C, pE, prec_cap,
                          (verbose - 1 if verbose > 0 else verbose))
    apf = [nfs[i].trace_of_frobenius(pf[i], power=powf[i])
           for i in range(nE)]
    Lf = [nfs[i].coefficient_field() for i in range(nE)]
    Lcom = [common_embedding_field(LE[i], Lf[i])
            for i in range(nE)]
    LphiE = [LE[i].embeddings(Lcom[i])[0] for i in range(nE)]
    Lphif = [nfs[i].embedding(Lcom[i]) for i in range(nE)]
    Bnew = ZZ(p * lcm(gcd([(phiE[i](apE[i]) - phif[i](apf[i])).absolute_norm()
                           for i in range(nE)]) for apE in apE_ls))
    Bnew = gcd(Bold, Bnew)
    if B != 0:
        Bnew = lcm(Bnew, ZZ(Bold / product(prime^(Bold.ord(prime))
                                           for prime in B)))
    return Bnew

def _eliminate_by_trace(curves, newforms, prime, B, C, prec_cap,
                        verbose):
    """An implementation of :func:`eliminate_by_trace`

    INPUT:

    - ``curves`` -- The Frey curves

    - ``newforms`` -- The newform list in the final format

    - ``prime`` -- The base prime over which to do elimination

    - ``B`` -- Zero or a list of prime numbers, the primes to be
      eliminated.

    - ``C`` -- The condition to be used

    - ``prec_cap`` -- The precision cap on the variables

    - ``verbose`` -- Verbosity argument

    """
    p = (prime if prime in ZZ and ZZ(prime).is_prime()
         else prime.smallest_integer())
    if len(newforms) == 0:
        return newforms
    if isinstance(newforms, ConditionalValue):
        return apply_to_conditional_value(lambda nfs, con:
                                          _eliminate_by_trace(curves, nfs, p,
                                                              B, con & C,
                                                              prec_cap,
                                                              verbose),
                                          newforms, use_condition=True,
                                          default_condition=C)
    if verbose > 0:
        print("Comparing traces of frobenius at " + str(p) + " for " +
               str(len(newforms)) + " cases.")
    nE = len(curves)
    KE = tuple((QQ if isinstance(curve, Qcurve) else curve.definition_field())
               for curve in curves)
    LE = tuple((curve.splitting_image_field() if isinstance(curve, Qcurve)
                else QQ) for curve in curves)
    pE = tuple((p if K == QQ else K.prime_above(p)) for K in KE)
    result = []
    Bprod = (B if B in ZZ else product(B))
    for nfs in newforms:
        Bnew = _single_elimination(curves, KE, LE, nfs, p, prime, pE,
                                   B, Bprod, C, prec_cap, verbose)
        if abs(Bnew) != 1:
            result.append(tuple([nfs[i] for i in range(nE)] + [Bnew]))
    _init_traces.clear_cache()
    return result

def eliminate_by_traces(curves, newforms, condition=None, primes=50,
                        precision_cap=1, verbose=False):
    r"""Eliminate newforms associated to Frey curves by the traces of
    Frobenius.

    Let $E$ be a Frey curve, $f$ be a newform of which the mod $l$
    galois representations might agree with the mod $l$ galois
    representation of $E$ and $F$ be a frobenius element of the common
    domain of the afore mentioned representations. By comparing the
    traces of these representations at $F$ we can check if the
    representations are indeed the same. If not we eliminate the
    newform for mod $l$ representations. This function will perform
    this elimination for multiple Frobenius elements.

    Note that the value of the trace of Frobenius for the Frey curve
    might depend on the value of the variables. In this case we
    eliminate the newform if the trace of Frobenius for the newform
    can not correspond to any of the possible traces for the elliptic
    curve. If we have multiple Frey curves we compare traces
    simultaneously for each combination of associated newforms, in
    which case we eliminate a tuple of newforms if there is no
    associated possible tuple of traces.

    INPUT:

    - ``curves`` -- A Frey curve or a list/tuple of Frey curves that
      share the same parameters.

    - ``newforms`` -- A list of newforms that are associated to the
      given frey curve, i.e. such that their mod-l representations
      could be isomorphic to the mod-l representation of the frey
      curve.  Instead of a list one can also give the value None, in
      which case the list will be retrieved from the curve by the
      method :meth:`newform_candidates` of the given Frey curve.

      If multiple Frey curves are given this argument should be a
      tuple of length equal to the number of given Frey curves,
      containing for each Frey curve an entry as described above in
      the same order as the Frey curves are given.

      Alternatively the input may also be a list of tuples wherein
      each tuple has length equal to the number of Frey curves and
      each entry thereof is a newform associated to the corresponding
      Frey curve.

      In the last input format, each tuple may also be 1 longer than
      the number of Frey curves, in which case the last entry must be
      an integer divisible by all primes l for which the mod-l
      representations of the given elliptic curves may still be
      isomorphic to the corresponding newforms in this tuple. This can
      be used to chain multiple elimination methods.

      Everywhere where a list is expected, one may also give a
      ConditionalValue of which the possible values are lists of the
      expected format. The method will be applied to each entry
      individually and the condition will be replaced by the condition
      that both the given condition and the condition at which that
      value is attained hold.

    - ``condition`` -- A Condition giving the restrictions on the
      parameters of the given Frey curve(s) that should be
      considered. By default this will be set to the condition stored
      in the first given Frey curve.

    - ``primes`` -- A list of finite primes of the field in which the
      parameters of the Frey curves live or a strictly positive
      integer (default: 50). The latter corresponds to all the finite
      primes of that field above a prime number smaller than the given
      integer. These will be the primes underlying the Frobenius
      elements for which the traces of the galois representations at
      that element should be compared. Each prime should be give as a
      prime number if the field of the parameters is the rationals and
      a prime ideal otherwise. Note that the Frobenius elements might
      be chosen for primes lying above these primes if the base fields
      of the Galois representations considered are larger than the
      field of the parameters.

    - ``precision_cap`` -- A non-negative integer (default: 1) giving
      the maximal precision used to compute the reduction type of a
      Frey curve at a prime. Note that setting this value to a value
      higher than 1 might result in a long computation time and heavy
      memory usage and is strongly discouraged.

    - ``verbose`` -- A boolean value or an integer (default:
      False). When set to True or any value larger then zero will
      print comments to stdout about the computations being done
      whilst busy. If set to False or 0 will not print such comments.
      If set to any negative value will also prevent the printing of
      any warnings.  A higher value will cause more messages to be
      printed.

    OUTPUT:
    
    A list of tuples, wherein each tuple contains:

    - for each frey curve given a corresponding newform for which an
      isomorphism between the mod $l$ galois representations of that
      curve and that newform could not be excluded by this method.

    - as its last entry an integer that is divisible by all prime
      numbers $l$ for which the mentioned mod $l$ galois
      representations could still be isomorphic. Note that this number
      is always 0 if there is infinitely many such $l$ remaining.

    .. SEEALSO::

        :func:`eliminate_by_trace`,
        :func:`kraus_method`,
        :func:`eliminate_cm_forms`,
        :func:`eliminate_primes`

    """
    curves, newforms, condition = _init_elimination_data(curves, newforms,
                                                         condition)
    K = curves[0]._R.fraction_field()
    if primes in ZZ and primes > 0:
        primes = prime_range(primes)
        if K != QQ:
            primes = [P for p in primes for P in K.primes_above(p)]
    if not (isinstance(primes, list) and
            all((K == QQ and p in ZZ and p > 0 and p.is_prime()) or
                (is_Ideal(p) and p in K.ideal_monoid() and p.is_prime())
                for p in primes)):
        raise ValueError("%s is not a list of valid primes"%(primes,))
    return apply_to_conditional_value(lambda nfs, con:
                                      _eliminate_by_traces(curves, nfs,
                                                           con & condition,
                                                           primes,
                                                           precision_cap,
                                                           verbose),
                                      newforms, use_condition=True,
                                      default_condition=condition)
        
def _eliminate_by_traces(curves, newforms, condition, primes, precision_cap,
                         verbose):
    r"""Implementation of :func:`eliminate_by_traces`.

    """
    for prime in primes:
        newforms = _eliminate_by_trace(curves, newforms, prime, 0, condition,
                                       precision_cap, verbose)
    return newforms

def kraus_method(curves, newforms, l, polynomials, primes=200, condition=None,
                 precision_cap=1, verbose=False):
    r"""Eliminate newforms associated to Frey curves by Kraus's method.

    Let $E$ be a Frey curve, $f$ be a newform of which the mod $l$
    galois representations might agree with the mod $l$ galois
    representation of $E$ and $F$ be a frobenius element of the common
    domain of the afore mentioned representations. By comparing the
    traces of these representations at $F$ we can check if the
    representations are indeed the same. If not we eliminate the
    newform for mod $l$ representations.

    Note that the value of the trace of Frobenius for the Frey curve
    might depend on the value of the variables. In this case we
    eliminate the newform if the trace of Frobenius for the newform
    can not correspond to any of the possible traces for the elliptic
    curve. If we have multiple Frey curves we compare traces
    simultaneously for each combination of associated newforms, in
    which case we eliminate a tuple of newforms if there is no
    associated possible tuple of traces.

    If some polynomial in the common variables of the Frey curve(s) is
    equal to an $l$-th power for $l$ prime, then this gives
    limitations for the possible residues of the variables in a
    residue field which has a number of elements congruent to $1$
    modulo $l$. For a fixed prime number $l$ this function computes
    the possible primes for which this happens and compares the traces
    of the associated Frobenius elements to eliminate newforms.

    INPUT:

    - ``curves`` -- A Frey curve or a list/tuple of Frey curves that
      share the same parameters.

    - ``newforms`` -- A list of newforms that are associated to the
      given Frey curve, i.e. such that their mod-l representations
      could be isomorphic to the mod-l representation of the Frey
      curve.  Instead of a list one can also give the value None, in
      which case the list will be retrieved from the curve by the
      method :meth:`newform_candidates` of the given Frey curve.

      If multiple Frey curves are given this argument should be a
      tuple of length equal to the number of given Frey curves,
      containing for each Frey curve an entry as described above in
      the same order as the Frey curves are given.

      Alternatively the input may also be a list of tuples wherein
      each tuple has length equal to the number of Frey curves and
      each entry thereof is a newform associated to the corresponding
      Frey curve.

      In the last input format, each tuple may also be 1 longer than
      the number of Frey curves, in which case the last entry must be
      an integer divisible by all primes l for which the mod-l
      representations of the given elliptic curves may still be
      isomorphic to the corresponding newforms in this tuple. This can
      be used to chain multiple elimination methods.

      Everywhere where a list is expected, one may also give a
      ConditionalValue of which the possible values are lists of the
      expected format. The method will be applied to each entry
      individually and the condition will be replaced by the condition
      that both the given condition and the condition at which that
      value is attained hold.

    - ``l`` -- The prime number $l$.

    - ``polynomials`` -- A polynomial or list thereof in the common
      variables of the given Frey curves, such that these polynomials
      are l-th powers in all possible values of the variables.

    - ``primes`` -- A list of finite primes of the field in which the
      parameters of the Frey curves live or a strictly positive
      integer (default: 200). The latter corresponds to all the finite
      primes of that field above a prime number smaller than the given
      integer. These will be the primes underlying the Frobenius
      elements for which the traces of the galois representations at
      that element should be compared. Each prime should be give as a
      prime number if the field of the parameters is the rationals and
      a prime ideal otherwise. Note that the Frobenius elements might
      be chosen for primes lying above these primes if the base fields
      of the Galois representations considered are larger than the
      field of the parameters. Note that only those primes are
      considered for which the associated residue field has order 1
      modulo $l$.

    - ``condition`` -- A Condition giving the restrictions on the
      variables of the given Frey curve(s) that should be
      considered. By default this will be set to the condition stored
      in the first given Frey curve.

    - ``precision_cap`` -- A non-negative integer (default: 1) giving
      the maximal precision used to compute the reduction type of a
      Frey curve at a prime. Note that setting this value to a value
      higher than 1 might result in a long computation time and heavy
      memory usage and is strongly discouraged.

    - ``verbose`` -- A boolean value or an integer (default:
      False). When set to True or any value larger then zero will
      print comments to stdout about the computations being done
      whilst busy. If set to False or 0 will not print such comments.
      If set to any negative value will also prevent the printing of
      any warnings.  A higher value will cause more messages to be
      printed.

    OUTPUT:
    
    A list of tuples, wherein each tuple contains:

    - for each frey curve given a corresponding newform for which an
      isomorphism between the mod $l$ galois representations of that
      curve and that newform could not be excluded by this method.

    - as its last entry an integer that is divisible by all prime
      numbers $l$ for which the mentioned mod $l$ galois
      representations could still be isomorphic. Note that this number
      is always 0 if there is infinitely many such $l$ remaining.

    .. SEEALSO::

        :func:`eliminate_by_trace`,
        :func:`eliminate_by_traces`,
        :func:`eliminate_cm_forms`,
        :func:`eliminate_primes`

    """
    curves, newforms, condition = _init_elimination_data(curves, newforms,
                                                         condition)
    if not (l in ZZ and l > 0 and l.is_prime()):
        raise ValueError("%s is not a prime number"%(prime,))
    if not isinstance(polynomials, list):
        if isinstance(polynomials, tuple):
            polynomials = list(polynomials)
        else:
            polynomials = [polynomials]
    K = curves[0]._R.fraction_field()
    if primes in ZZ and primes > 0:
        primes = prime_range(primes)
        if K != QQ:
            primes = [P for p in primes for P in K.primes_above(p)]
    if not (isinstance(primes, list) and
            all((K == QQ and p in ZZ and p > 0 and p.is_prime()) or
                (is_Ideal(p) and p in K.ideal_monoid() and p.is_prime())
                for p in primes)):
        raise ValueError("%s is not a list of valid primes"%(primes,))
    return apply_to_conditional_value(lambda nfs, con:
                                      _kraus_method(curves, nfs, l,
                                                    polynomials, primes,
                                                    con & condition,
                                                    precision_cap, verbose),
                                      newforms, use_condition=True,
                                      default_condition=condition)

def _kraus_method(curves, newforms, l, polynomials, primes, condition,
                  precision_cap, verbose):
    r"""Implementation of :func:`kraus_method`.

    """
    fields = tuple(poly.base_ring() for poly in polynomials)
    for p in primes:
        Ps = [_primes_1_mod_l(K, p, l) for K in fields]
        for P in itertools.product(*Ps):
            CP = [_init_kraus_condition(polynomials[i], P[i], l)
                  for i in range(len(P))]
            C = reduce(lambda x, y: x & y, CP, condition)
            newforms = _eliminate_by_trace(curves, newforms, p, [l], C,
                                           precision_cap, verbose)
    return newforms

def _primes_1_mod_l(K, p, l):
    r"""Find the primes in $K$ above $p$ that are 1 mod $l$.

    Find the primes $P$ in a number field $K$ above a prime number $p$
    such that the residue field of $P$ has a number of elements that
    is $1$ modulo $l$

    INPUT:

    - ``K`` -- The number field $K$

    - ``p`` -- The prime number $p$

    - ``l`` -- The prime number $l$

    OUTPUT:

    A list of primes of $K$ above $p$ for which the number of elements
    of the residue field is $1$ modulo $l$.

    """
    if K == QQ:
        if mod(p, l) == 1:
            return [p]
        else:
            return []
    else:
        return [P for P in K.primes_above(p)
                if mod(P.residue_field().cardinality(), l) == 1]
                        
def _init_kraus_condition(polynomial, prime, l):
    r"""Get the condition on the variable given by Kraus's method.

    INPUT:

    - ``polynomial`` -- A polynomial whose possible values are $l$-th
      powers.

    - ``prime`` -- A finite prime of a number field.

    - ``l`` -- The prime number $l$

    OUTPUT:

    The condition on the variables given by the fact that the values
    of the polynomial in the residue field of the given prime should
    be $l$-th powers.

    """
    if prime in ZZ:
        F = Integers(prime)
        u = F.unit_gens()[0]
    else:
        F = prime.residue_field()
        u = F.primitive_element()
    values = [F.lift(0)] + [F.lift(u^(n*l))
                            for n in range((F.cardinality() - 1)/l)]
    conditions = [CongruenceCondition(polynomial - val, prime)
                  for val in values]
    def combine_conditions(C1, C2):
        if C1 is None:
            return C2
        elif C2 is None:
            return C1
        else:
            return C1 | C2
    return reduce(combine_conditions, conditions, None)

def eliminate_cm_forms(curves, newforms, has_cm=False, condition=None):
    r"""Eliminate CM newforms associated to non-CM Frey curves.

    This function eliminates all those tuples of newforms for which
    one of the newforms has complex multiplication, but the
    corresponding Frey curve does not.

    INPUT:

    - ``curves`` -- A Frey curve or a list/tuple of Frey curves that
      share the same parameters.

    - ``newforms`` -- A list of newforms that are associated to the
      given frey curve, i.e. such that their mod-l representations
      could be isomorphic to the mod-l representation of the frey
      curve.  Instead of a list one can also give the value None, in
      which case the list will be retrieved from the curve by the
      method :meth:`newform_candidates` of the given Frey curve.

      If multiple Frey curves are given this argument should be a
      tuple of length equal to the number of given Frey curves,
      containing for each Frey curve an entry as described above in
      the same order as the Frey curves are given.

      Alternatively the input may also be a list of tuples wherein
      each tuple has length equal to the number of Frey curves and
      each entry thereof is a newform associated to the corresponding
      Frey curve.

      In the last input format, each tuple may also be 1 longer than
      the number of Frey curves, in which case the last entry must be
      an integer divisible by all primes l for which the mod-l
      representations of the given elliptic curves may still be
      isomorphic to the corresponding newforms in this tuple. This can
      be used to chain multiple elimination methods.

      Everywhere where a list is expected, one may also give a
      ConditionalValue of which the possible values are lists of the
      expected format. The method will be applied to each entry
      individually and the condition will be replaced by the condition
      that both the given condition and the condition at which that
      value is attained hold.

    - ``has_cm`` -- A boolean value or a list thereof (default:
      False). If it is a list it must have length equal to the number
      of considered Frey curves. Each entry of this list should be
      False if it is certain that the corresponding Frey curve does
      not have complex multiplication and True otherwise. If the given
      argument is a boolean value it will be turned into a list of the
      required length with each entry equal to the given value.

    - ``condition`` -- A Condition giving the restrictions on the
      parameters of the given Frey curve that should be considered. By
      default this will be set to the condition stored in the first
      given Frey curve. Note that this is only relevant if the
      newforms are not given yet.

    OUTPUT:
    
    A list of tuples, wherein each tuple contains:

    - for each frey curve given a corresponding newform for which an
      isomorphism between the mod $l$ galois representations of that
      curve and that newform could not be excluded by this method.

    - as its last entry an integer that is divisible by all prime
      numbers $l$ for which the mentioned mod $l$ galois
      representations could still be isomorphic. Note that this number
      is always 0 if there is infinitely many such $l$ remaining.

    .. SEEALSO::

        :func:`eliminate_by_trace`,
        :func:`eliminate_by_traces`,
        :func:`kraus_method`,
        :func:`eliminate_primes`

    """
    curves, newforms, condition = _init_elimination_data(curves, newforms,
                                                         condition)
    if isinstance(has_cm, tuple):
        has_cm = list(has_cm)
    if not isinstance(has_cm, list):
        has_cm = [has_cm for i in range(len(curves))]
    if len(has_cm) != len(curves):
        raise ValueError("Expected cm information for " + str(len(curves)) +
                         " curves, but got it for " + str(len(has_cm)) +
                         " curves")
    return _eliminate_cm_forms(newforms, has_cm)

def _eliminate_cm_forms(newforms, has_cm):
    r"""The implementation of :func:`eliminate_cm_forms`.

    """
    if isinstance(newforms, ConditionalValue):
        return apply_to_conditional_value(lambda nfs:
                                          _eliminate_cm_forms(nfs, has_cm),
                                          newforms)
    result = []
    for nfs in newforms:
        keep = True
        for i in range(len(has_cm)):
            keep = has_cm[i] or not nfs[i].has_cm()
            if not keep:
                break
        if keep:
            result.append(nfs)
    return result

def eliminate_primes(curves, newforms, N, condition=None):
    r"""Eliminate newforms associated to Frey curves for certain primes.

    For multiple prime numbers $l$ this function eliminates all those
    newforms for which the mod $l$ galois representation might agree
    with the mod $l$ galois representation of the associated Frey
    curve.

    INPUT:

    - ``curves`` -- A Frey curve or a list/tuple of Frey curves that
      share the same parameters.

    - ``newforms`` -- A list of newforms that are associated to the
      given frey curve, i.e. such that their mod-l representations
      could be isomorphic to the mod-l representation of the frey
      curve.  Instead of a list one can also give the value None, in
      which case the list will be retrieved from the curve by the
      method :meth:`newform_candidates` of the given Frey curve.

      If multiple Frey curves are given this argument should be a
      tuple of length equal to the number of given Frey curves,
      containing for each Frey curve an entry as described above in
      the same order as the Frey curves are given.

      Alternatively the input may also be a list of tuples wherein
      each tuple has length equal to the number of Frey curves and
      each entry thereof is a newform associated to the corresponding
      Frey curve.

      In the last input format, each tuple may also be 1 longer than
      the number of Frey curves, in which case the last entry must be
      an integer divisible by all primes l for which the mod-l
      representations of the given elliptic curves may still be
      isomorphic to the corresponding newforms in this tuple. This can
      be used to chain multiple elimination methods.

      Everywhere where a list is expected, one may also give a
      ConditionalValue of which the possible values are lists of the
      expected format. The method will be applied to each entry
      individually and the condition will be replaced by the condition
      that both the given condition and the condition at which that
      value is attained hold.

    - ``N`` -- A non-zero integer, divisible by all those primes $l$
      for which the newforms that could still be related by mod $l$
      galois representations should be eliminated.

    - ``condition`` -- A Condition giving the restrictions on the
      parameters of the given Frey curve(s) that should be
      considered. By default this will be set to the condition stored
      in the first given Frey curve. Note that this is only relevant
      if the newforms are not given yet.

    OUTPUT:
    
    A list of tuples, wherein each tuple contains:

    - for each frey curve given a corresponding newform for which an
      isomorphism between the mod $l$ galois representations of that
      curve and that newform could not be excluded by this method.

    - as its last entry an integer that is divisible by all prime
      numbers $l$ for which the mentioned mod $l$ galois
      representations could still be isomorphic. Note that this number
      is always 0 if there is infinitely many such $l$ remaining.

    .. SEEALSO::

        :func:`eliminate_by_trace`,
        :func:`eliminate_by_traces`,
        :func:`kraus_method`,
        :func:`eliminate_cm_forms`

    """
    curves, newforms, condition = _init_elimination_data(curves, newforms,
                                                         condition)
    if N not in ZZ or N == 0:
        raise ValueError("Argument N is not a non-zero integer.")
    return _eliminate_primes(newforms, N.prime_factors())

def _eliminate_primes(newforms, ls):
    r"""Implementation of :func:`eliminate_primes`.

    """
    if isinstance(newforms, ConditionalValue):
        return apply_to_conditional_value(lambda nfs:
                                          _eliminate_primes(nfs, ls), newforms)
    result = []
    for nfs in newforms:
        if nfs[-1] != 0:
            nfs = list(nfs)
            for l in ls:
                nfs[-1] = ZZ(nfs[-1] / (l^nfs[-1].ord(l)))
            nfs = tuple(nfs)
        if abs(nfs[-1]) != 1 :
            result.append(nfs)
    return result

def combine_newforms(*newforms):
    r"""Combine the output of different newform elimination methods into
    one.
    
    Given multiple lists of newforms given as output by one of the
    functions :func:`eliminate_by_traces`, :func:`eliminate_by_trace`,
    :func:`kraus_method`, :func:`eliminate_cm_forms` or
    :func:`eliminate_primes` outputs a single list of tuples of
    newforms that combines all this information.

    INPUT:
    
    Any number of arguments, wherein each argument is a list of tuples
    of the form $(f_1, ..., f_m, N)$, where each $f_i$ is a newform
    and $N$ is an integer. The integer $m$ should be the same among
    all tuples in a list.

    OUTPUT:
    
    A list of tuples, such that for each combination of tuples $(f_(1,
    1), ..., f_(1, m_1), N_1), ..., (f_(n, 1), ..., f_(n, m_n), N_n)$
    from the respective lists given, the tuple $(f_(1, 1), ..., f(1,
    m_1), f_(2, 1), ..., f_(n, m_n), gcd(N_1, ..., N_n))$ is in this
    list. Note that tuples wherein the last entry is 1 or -1 are
    omitted.

    """
    newforms = conditional_product(*newforms)
    return apply_to_conditional_value(lambda nfs: _combine_newforms(*nfs),
                                      newforms)

def _combine_newforms(*newforms):
    r"""Implementation of :func:`combine_newforms`"""
    return [nfs for nfs in (sum((item[:-1] for item in items), ()) + 
                            (gcd([item[-1] for item in items]),)
                            for items in itertools.product(*newforms))
            if abs(nfs[-1]) != 1]
