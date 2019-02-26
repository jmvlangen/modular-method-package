r"""An implementation of Tate's algorithm for Frey curves

This code allows for the computation of the possible conductors of an
elliptic curve which depends on some integral parameters.

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

from sage.structure.sage_object import SageObject

from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class

from sage.all import Infinity

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
        
def _least_power(poly_list, pAdics):
    result = Infinity
    for poly in poly_list:
        result = min([result]+[pAdics.valuation(c)
                               for c in poly.coefficients()])
    return result

def _determine_level(poly_list, pAdics, T, max_value):
    level = max_value - _least_power(poly_list, pAdics)
    if level > -Infinity:
        level = ceil(level/pAdics.extension_multiplicity(T.pAdics()))
    if level < 0:
        level = 0
    return level
        
def _number_of_zero_entries(iterable):
    count = 0
    for c in iterable:
        if c == 0:
            count += 1
    return count

def get_cases_invariant(poly, pAdics, T, name, general_case, variables=None,
                        verbose=False, precision=20, result=[],
                        precision_cap=20, **kwds):
    Tlist, v_min = find_pAdic_roots(poly, pAdics=pAdics,
                                    variables=variables, value_tree=T,
                                    precision=20, verbose=verbose,
                                    give_list=True,
                                    precision_cap=precision_cap)
    for i in range(len(Tlist)):
        if not Tlist[i].is_empty():
            case = general_case.copy()
            case[name] = v_min + i
            case['T'] = Tlist[i]
            result.append(case)
    return result
        

def get_two_cases_invariant(poly, pAdics, T, bndry, case_big, case_small,
                            variables=None, verbose=False, result=[],
                            precision_cap=20, **kwds):
    Tbig, Tsmall = find_pAdic_roots(poly, pAdics=pAdics,
                                    variables=variables, value_tree=T,
                                    precision=bndry, verbose=verbose,
                                    precision_cap=precision_cap)
    if not Tbig.is_empty():
        case_big['T'] = Tbig
        result.append(case_big)
    if not Tsmall.is_empty():
        case_small['T'] = Tsmall
        result.append(case_small)
    return result
        
def get_number_of_roots_cases(poly_list, pAdics, T, name, general_case,
                              variables=None, verbose=False, result=[],
                              **kwds):
    r"""
    Let n be the length of poly_list. For each tuple in T this will
    determine how many roots the polynomial
        poly_list[0]*X^n + ... + poly_list[n-2]*X + poly_list[n-1]
    has over the residue field corresponding to T. It will store this
    information in a copy of case under the given name.
    """
    Tdict = {}
    n = _determine_level(poly_list, pAdics, T, 1)
    F = pAdics.residue_field()
    S.<x> = F[]
    child_list = T.children_at_level(n)
    if verbose > 0:
        print "Checking irreducibility in %d cases"%len(child_list)
    for child in child_list:
        coeff_list = [F(poly(child.representative())) for poly in poly_list]
        f = 0
        k = len(coeff_list)-1
        for i in range(k+1):
            f += coeff_list[i] * x^(k-i)
        m = sum([r[1] for r in f.roots()])
        if not Tdict.has_key(m):
            Tdict[m] = pAdicNode(pAdics=T.pAdics(), width=T.width)
        Tdict[m].merge(child, from_root=True)
    for (m,Tm) in Tdict.iteritems():
        case = general_case.copy()
        case[name] = m
        case['T'] = Tm
        result.append(case)
    return result
    
def _tate_step_1(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=2, E0=E0)
    case_small = dict(E=E, vDelta=0, m=1, f=0, c=1, KS="I0", E0=E0) 
    return get_two_cases_invariant(S(E.discriminant()), pAdics, T, 1,
                                   case_big, case_small, **kwds)
    
def _tate_step_2_transformation(E, S, pAdics, T, E0, verbose=False, result=[],
                                **kwds):
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    char = pAdics.characteristic()
    replaceCases = dict()
    
    # First transforming the singular point to (0,0), which (could) give us different cases.
    s = 0
    if char == 2:
        if R == ZZ and P in ZZ.ideal_monoid():
            s = ((ZZ.quotient(P.gen() - 1)(2))^(-1)).lift()
        else:
            s = ((ZZ.quotient(ZZ(norm(P) - 1))(2))^(-1)).lift()
    if char == 3:
        if R == ZZ and P in ZZ.ideal_monoid():
            s = ((ZZ.quotient(P.gen() - 1)(3))^(-1)).lift() 
        else:
            s = ((ZZ.quotient(ZZ(norm(P) - 1))(3))^(-1)).lift()
    if s == 0:
        s = 1
    
    level = _determine_level([S(E.a1()), S(E.a2()), S(E.a3()), S(E.a4()),
                              S(E.a6())], pAdics, T, 1)
    if verbose > 0:
        print "Determining singular point for %d cases"%(T.count_children_at_level(level),)
    for node in T.children_at_level(level):
        # Isolating the singular point x,y
        a1 = F(S(E.a1())(node.representative()))
        a2 = F(S(E.a2())(node.representative()))
        a3 = F(S(E.a3())(node.representative()))
        a4 = F(S(E.a4())(node.representative()))
        a6 = F(S(E.a6())(node.representative()))
    
        if char == 2:
            if a1 == 0:
                x = a4^s
                y = (a6 + a2 * a4)^s
            else:
                x = a3 / a1
                y = (x^2 + a4) / a1
        else:
            # Coordinate transformation to end up in the case where a1 = a3 = 0,
            # hence y = 0
            a22 = a2 + a1^2 / 4
            a42 = a4 + a1 * a3 / 2
            a62 = a6 + a3^2 / 4
            y = 0
            
            if char == 3:
                if a22 == 0:
                    x = (-a62)^s
                else:
                    x = - a42 / a22
            else:
                # Coordinate transformation to end up in the case where also a2 = 0
                a43 = a42 - a22^2 / 3
                a63 = a62 - a22 * a42 / 3 + 2 * a22^3 / 27
                
                if a43 == 0:
                    x = 0
                else:
                    x = - 3 * a63 / ( 2 * a43 )
                
                # Transforming back
                x = x - a22 / 3
                
            # Transforming back
            y = y - a1 * x / 2 - a3 / 2
        
        singularPoint = tuple([x,y])
        if singularPoint in replaceCases:
            replaceCases[singularPoint].merge(node, from_root=True)
        else:
            Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
            Tn.merge(node, from_root=True)
            replaceCases[singularPoint] = Tn
    
    if verbose > 0:
        print "Performing transformation for %d cases"%(len(replaceCases),)
    for (point,Tn) in replaceCases.iteritems():
        xn = F.lift(point[0])
        yn = F.lift(point[1])
        En = E.rst_transform(xn,0,yn)
        result.append(dict(next_step=2+1/2, T=Tn, E=En, E0=E0))
            
    return result
    
def _tate_step_2(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=3, E0=E0)
    case_small = dict(E=E, f=1, KS="In", E0=E0) 
    return get_two_cases_invariant(S(E.b2()), pAdics, T, 1,
                                   case_big, case_small, **kwds)

def _tate_step_3(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=4, E0=E0)
    case_small = dict(E=E, KS="II", m=1, c=1, E0=E0)
    char = pAdics.characteristic()
    if char != 2 and char != 3:
        case_small['vDelta'] = 2
        case_small['f'] = 2
    return get_two_cases_invariant(S(E.a6()), pAdics, T, 2,
                                   case_big, case_small, **kwds)
    
def _tate_step_4(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=5, E0=E0)
    case_small = dict(E=E, KS="III", m=2, c=2, E0=E0)
    char = pAdics.characteristic()
    if char != 2:
        case_small['vDelta'] = 3
        case_small['f'] = 2
    return get_two_cases_invariant(S(E.b8()), pAdics, T, 3,
                                   case_big, case_small, **kwds)

def _tate_step_5(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=6, E0=E0)
    case_small = dict(E=E, KS="IV", m=3, E0=E0)
    char = pAdics.characteristic()
    if char != 3:
        case_small['vDelta'] = 4
        case_small['f'] = 2
    return get_two_cases_invariant(S(E.b6()), pAdics, T, 3,
                                   case_big, case_small, **kwds)

def _tate_step_6_transformation(E, S, pAdics, T, E0, verbose=False, result=[],
                                **kwds):
    pi = pAdics.uniformizer()
    char = pAdics.characteristic()
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    
    changeDict = dict()   
    if char != 2:
        a1 = S(E.a1())
        a31 = S(E.a3()/pi)     
        half = (F(2)^(-1))
        
        level = _determine_level([a1, a31], pAdics, T, 1)
        if verbose > 0:
            print "Determining necessary transformation for %d cases"%(T.count_children_at_level(level),)
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
            sqrtPower = ((ZZ.quotient(P.gen() - 1)(2))^(-1)).lift()
        else:
            sqrtPower = ((ZZ.quotient(ZZ(norm(P) - 1))(2))^(-1)).lift()
        if sqrtPower == 0:
            sqrtPower = 1
            
        a2 = S(E.a2())
        a62 = S(E.a6()/(pi^2))
        
        level = _determine_level([a2, a62], pAdics, T, 1)
        if verbose > 0:
            print "Determining necessary transformation for %d cases"%T.count_children_at_level(level)
        for node in T.children_at_level(level):
            alphaSqrd = -F(a2(node.representative()))
            betaSqrd = F(-a62(node.representative()))
            alpha = alphaSqrd^sqrtPower
            beta = betaSqrd^sqrtPower
            alphaBetaPair = tuple([alpha, beta])
            if alphaBetaPair in changeDict:
                changeDict[alphaBetaPair].merge(node, from_root=True)
            else:
                Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
                Tn.merge(node, from_root=True)
                changeDict[alphaBetaPair] = Tn
                
    if verbose > 0:
        print "Performing %d transformations"%len(changeDict)
    for (alphaBetaPair, Tn) in changeDict.iteritems():
        En = E.rst_transform(0, F.lift(alphaBetaPair[0]),F.lift(alphaBetaPair[1]) * pi)
        result.append(dict(next_step=6+1/2, T=Tn, E=En, E0=E0))
    return result

def _tate_step_6(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=7, E0=E0)
    case_small = dict(E=E, KS="I0*", m=5, E0=E0)
    char = pAdics.characteristic()
    if char != 2:
        case_small['vDelta'] = 6
        case_small['f'] = 2
    return get_two_cases_invariant(S(-4*E.a2()^3*E.a6() + E.a2()^2*E.a4()^2 - 4*E.a4()^3 - 27*E.a6()^2 + 18*E.a2()*E.a4()*E.a6()),
                                   pAdics, T, 7, case_big, case_small, **kwds)
        
def _tate_step_7(E, S, pAdics, T, E0, case, restrictions, **kwds):
    case_big = dict(E=E, next_step=8, E0=E0)
    if pAdics.characteristic() == 2 or \
       'type' in restrictions or \
       _should_calculate_c(case, restrictions):
        case_small = dict(E=E, next_step=7+encodeQuotient(1,0), E0=E0)
    else:
        case_small = dict(E=E, KS="In*", f=2, E0=E0)
    return get_two_cases_invariant(S(3*E.a4() - (E.a2())^2), pAdics, T, 3,
                                   case_big, case_small, **kwds)

def encodeQuotient(n, b):
    m = 2*n + b
    return 1 / (m+1)

def decodeQuotient(q):
    m = ZZ(1/q) - 1
    if is_even(m):
        n = ZZ(m/2)
        b = 0
    else:
        n = ZZ((m - 1)/2)
        b = 1   
    return n , b
            
def _tate_step_7sub_transformation(E, S, pAdics, T, E0, n, verbose=False,
                                   result=[], **kwds):
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    pi = pAdics.uniformizer()
    char = pAdics.characteristic()
    
    if char == 2:
        if R == ZZ and P in ZZ.ideal_monoid():
            sqrtPower = ((ZZ.quotient(P.gen() - 1)(2))^(-1)).lift()
        else:
            sqrtPower = ((ZZ.quotient(ZZ(norm(P) - 1))(2))^(-1)).lift()
        if sqrtPower == 0:
            sqrtPower = 1
    
    changeDict = dict()
    if n == 1:
        a42 = S(E.a4()/(pi^2))
        if char == 2:
            level = _determine_level([a42], pAdics, T, 1)
        else:
            a21 = S(E.a2()/pi)
            if char == 3:
                level = _determine_level([a21, a42], pAdics, T, 1)
            else:
                a63 = S(E.a6()/(pi^3))
                level = _determine_level([a21, a42, a63], pAdics, T, 1)
    else:
        if is_odd(n):
            k = ZZ((n+1)/2)
            a21 = S(E.a2()/pi)
            if char == 2:
                a6k = S(E.a6()/(pi^(2*k+1)))
                level = _determine_level([a21, a6k], pAdics, T, 1)
            else:
                a4k = S(E.a4()/(pi^(k+1)))
                level = _determine_level([a21, a4k], pAdics, T, 1)
        else:
            k = ZZ((n+2)/2)
            if char == 2:
                a6k = S(E.a6()/(pi^(2*k)))
                level = _determine_level([a6k], pAdics, T, 1)
            else:
                a3k = S(E.a3()/(pi^k))
                level = _determine_level([a3k], pAdics, T, 1)
      
    if verbose > 0:
        print "Determining necessary transformation for %d cases"%T.count_children_at_level(level)
    for node in T.children_at_level(level):
        if n == 1:
            a42_ev = F(a42(node.representative()))
            if char == 2:
                change = a42_ev^sqrtPower
            else:
                a21_ev = F(a21(node.representative()))
                if char == 3:
                    change = -a42_ev / (2 * a21_ev)
                else:
                    a63_ev = F(a63(node.representative()))
                    change = (a21_ev*a42_ev - 9*a63_ev)/(2*(3*a42_ev - a21_ev^2))
        else:
            if is_odd(n):
                a21_ev = F(a21(node.representative()))
                if char == 2:
                    a6k_ev = F(a6k(node.representative()))
                    square = a6k_ev / a21_ev
                    change = square^sqrtPower
                else:
                    a4k_ev = F(a4k(node.representative()))
                    change = - a4k_ev / (2 * a21_ev)
            else:
                if char == 2:
                    a6k_ev = F(a6k(node.representative()))
                    change = a6k_ev^sqrtPower
                else:
                    a3k_ev = F(a3k(node.representative()))
                    change = - a3k_ev / F(2)
                    
        if change in changeDict:
            changeDict[change].merge(node, from_root=True)
        else:
            Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
            Tn.merge(node, from_root=True)
            changeDict[change] = Tn
            
    if verbose > 0: 
        print "Performing %d transformations."%len(changeDict)
    for (change, Tn) in changeDict.iteritems():
        if n==1:
            En = E.rst_transform(pi * F.lift(change), 0, 0)
        elif is_odd(n):
            En = E.rst_transform(pi^k * F.lift(change), 0, 0)
        else:
            En = E.rst_transform(0, 0, pi^k * F.lift(change))
        result.append(dict(next_step=7 + encodeQuotient(n ,1), T=Tn, E=En, E0=E0))
        
    return result
        
def _tate_step_7sub(E, S, pAdics, T, E0, n, **kwds):
    case_big = dict(E=E, next_step=7+encodeQuotient(n+1,0), E0=E0)
    case_small = dict(E=E, KS="I"+str(n)+"*", m=5+n, E0=E0)
    if is_odd(n):
        k = ZZ((n+3)/2)
        return get_two_cases_invariant(S(E.a3()^2 + 4*E.a6()), pAdics, T,
                                       2*k+1, case_big, case_small, **kwds)
    else:
        k = ZZ((n+2)/2)
        return get_two_cases_invariant(S(E.a4()^2 - 4*E.a2()*E.a6()), pAdics, T,
                                       2*k + 3, case_big, case_small, **kwds)
     
def _tate_step_8_transformation(E, S, pAdics, T, E0, verbose=False, result=[],
                                **kwds):
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    pi = pAdics.uniformizer()
    
    if F.characteristic() != 3:
        result.append(dict(next_step=8+1/2, T=T,
                           E=E.rst_transform(-E.a2() * F.lift(F(3)^(-1)), 0, 0),
                           E0=E0))
    else:
        if R == ZZ and P in ZZ.ideal_monoid():
            cubertPower = ((ZZ.quotient(P.gen() - 1)(3))^(-1)).lift()
        else:
            cubertPower = ((ZZ.quotient(ZZ(norm(P) - 1))(3))^(-1)).lift()
        
        changeDict = dict()
        a63 = S(E.a6()/(pi^3))
        level = _determine_level([a63], pAdics, T, 1)
        if verbose > 0:
            print "Determining necessary transformation for %d cases"%T.count_children_at_level(level)
        for node in T.children_at_level(level):
            cube = F(a63(node.representative()))
            change = cube^cubertPower
            if change in changeDict:
                changeDict[change].merge(node, from_root=True)
            else:
                Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
                Tn.merge(node, from_root=True)
                changeDict[change] = Tn
        
        if verbose > 0:
            print "Performing %d transformations."%len(changeDict)
        for (change, Tn) in changeDict.iteritems():
            En = E.rst_transform(-pi * F.lift(change), 0, 0)
            result.append(dict(next_step=8+1/2, T=Tn, E=En, E0=E0))
        
    return result
            
def _tate_step_8(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=9, E0=E0)
    case_small = dict(E=E, KS="IV*", m=7, E0=E0)
    char = pAdics.characteristic()
    if char != 3:
        case_small['vDelta'] = 8
        case_small['f'] = 2
    return get_two_cases_invariant(S( ( E.a3() )^2 + 4 * E.a6() ), pAdics, T,
                                   5, case_big, case_small, **kwds)
        
def _tate_step_9_transformation(E, S, pAdics, T, E0, verbose=False, result=[],
                                **kwds):
    R = pAdics.order()
    P = pAdics.prime_ideal()
    F = pAdics.residue_field()
    pi = pAdics.uniformizer()
    
    if F.characteristic() != 2:
        result.append(dict(next_step=9+1/2, T=T,
                           E=E.rst_transform(0, 0, -E.a3() * F.lift(F(2)^(-1))), E0=E0))
    else:
        if R == ZZ and P in ZZ.ideal_monoid():
            sqrtPower = ( ( ZZ.quotient( P.gen() - 1 )(2) )^(-1) ).lift()
        else:
            sqrtPower = ( ( ZZ.quotient( ZZ(norm(P) - 1) )(2) )^(-1) ).lift()
        if sqrtPower == 0:
            sqrtPower = 1
        
        changeDict = dict()
        a64 = S(E.a6()/(pi^4))
        level = _determine_level([a64], pAdics, T, 1)
        if verbose > 0:
            print "Determining necessary transformation for %d cases"%T.count_children_at_level(level)
        for node in T.children_at_level(level):
            square = -F(a64(node.representative()))
            change = square^sqrtPower
            if change in changeDict:
                changeDict[change].merge(node, from_root=True)
            else:
                Tn = pAdicNode(pAdics=T.pAdics(), width=T.width)
                Tn.merge(node, from_root=True)
                changeDict[change] = Tn
        
        if verbose > 0:
            print "Performing %d transformations."%len(changeDict)
        for (change, Tn) in changeDict.iteritems():
            En = E.rst_transform(0, 0, -pi^2 * F.lift(change))
            result.append(dict(next_step=9+1/2, T=Tn, E=En, E0=E0))
        
    return result
            
def _tate_step_9(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=10, E0=E0)
    case_small = dict(E=E, KS="III*", m=8, c=2, E0=E0)
    char = pAdics.characteristic()
    if char != 2:
        case_small['vDelta'] = 9
        case_small['f'] = 2
    return get_two_cases_invariant(S( E.a4() ), pAdics, T, 4,
                                   case_big, case_small, **kwds)
            
def _tate_step_10(E, S, pAdics, T, E0, **kwds):
    case_big = dict(E=E, next_step=11, E0=E0)
    case_small = dict(E=E, KS="II*", m=9, c=1, E0=E0)
    char = pAdics.characteristic()
    if char != 2 and char != 3:
        case_small['vDelta'] = 10
        case_small['f'] = 2
    return get_two_cases_invariant(S(E.a6()), pAdics, T, 6,
                                   case_big, case_small, **kwds)
                               
def _tate_step_11(E, S, pAdics, T, E0, verbose=False, result=[], **kwds):
    if verbose > 0:
        print "Performing final transformation and restarting the algorithm."
    pi = pAdics.uniformizer()
    a1 = S(E.a1()/pi)
    a2 = S(E.a2()/(pi^2))
    a3 = S(E.a3()/(pi^3))
    a4 = S(E.a4()/(pi^4))
    a6 = S(E.a6()/(pi^6))
    E = EllipticCurve([a1,a2,a3,a4,a6])
    result.append(dict(next_step=1, T=T, E=E, E0=E))
    return result
    
def _tate_calculate_vDelta(E, S, pAdics, T, case, **kwds):
    return get_cases_invariant(S(E.discriminant()), pAdics, T, 'vDelta',
                               case, **kwds)
    
def _tate_calculate_m(E, S, pAdics, T, case, result=[], **kwds):
    KS = case['KS']
    if KS.endswith('*'):
        if KS.startswith('IV'): #IV*
            case['m'] = 7
        elif KS.startswith('III'): #III*
            case['m'] = 8
        elif KS.startswith('II'): #II*
            case['m'] = 9
        elif KS.startswith('I0'): #I0*
            case['m'] = 5
        else: #In*
            case['m'] = case['vDelta'] - 1
    else:
        if KS.startswith('IV'): #IV
            case['m'] = 3
        elif KS.startswith('III'): #III
            case['m'] = 2
        elif KS.startswith('II'): #II
            case['m'] = 1
        elif KS.startswith('I0'): #II
            case['m'] = 1
        else: #In
            case['m'] = case['vDelta']
    result.append(case)
    return result
    
def _tate_calculate_c(E, S, pAdics, T, case, result=[], **kwds):
    KS = case['KS']
    pi = pAdics.uniformizer()
    if KS.endswith('*'):
        if KS.startswith('IV'): #IV*
            if case.has_key('roots'):
                case['c'] = case['roots'] + 1
            else:
                return get_number_of_roots_cases([S(1),
                                                  S(E.a3()/(pi^2)),
                                                  S(-E.a6()/(pi^4))],
                                                 pAdics, T, 'roots', case,
                                                 result=result, **kwds)
        elif KS.startswith('III'): #III*
            case['c'] = 2
        elif KS.startswith('II'): #II*
            case['c'] = 1
        elif KS.startswith('I0'): #I0*
            if case.has_key('roots'):
                case['c'] = case['roots'] + 1
            else:
                return get_number_of_roots_cases([S(1),
                                                  S(E.a2()/pi),
                                                  S(E.a4()/(pi^2)),
                                                  S(E.a6()/(pi^3))],
                                                 pAdics, T, 'roots', case,
                                                 result=result, **kwds)
        else: #In*
            if case.has_key('roots'):
                case['c'] = case['roots'] + 2
            else:
                n = Integer(KS[1:-1])
                if is_odd(n):
                    k = ZZ((n+3)/2)
                    return get_number_of_roots_cases([S(1),
                                                      S(E.a3()/(pi^(k))),
                                                      S(-E.a6()/(pi^(2*k)))],
                                                     pAdics, T, 'roots', case,
                                                     result=result, **kwds)
                else:
                    k = ZZ((n+2)/2)
                    return get_number_of_roots_cases([S(E.a2()/pi),
                                                      S(E.a4()/(pi^(k+1))),
                                                      S(E.a6()/(pi^(2*k+1)))],
                                                     pAdics, T, 'roots', case,
                                                     result=result, **kwds)
    else:
        if KS.startswith('IV'): #IV
            if case.has_key('roots'):
                case['c'] = case['roots'] + 1
            else:
                return get_number_of_roots_cases([S(1),
                                                  S(E.a3()/pi),
                                                  S(-E.a6()/(pi^2))],
                                                 pAdics, T, 'roots', case,
                                                 result=result, **kwds)
        elif KS.startswith('III'): #III
            case['c'] = 2
        elif KS.startswith('II'): #II
            case['c'] = 1
        elif KS.startswith('I0'): #II
            case['c'] = 1
        else: #In
            if case['split']:
                case['c'] = case['vDelta']
            else:
                case['c'] = (1 if is_odd(case['vDelta']) else 2)

    result.append(case)
    return result
    
def _tate_calculate_f(E, S, pAdics, T, case, result=[], **kwds):
    case['f'] = case['vDelta'] - case['m'] + 1
    result.append(case)
    return result
    
def _tate_calculate_n(E, S, pAdics, T, case, result=[], **kwds):
    KS = case['KS']
    if KS.endswith('*'): #In*
        case['KS'] = KS.replace('n', str(case['vDelta']-6))
    else: #In
        case['KS'] = KS.replace('n', str(case['vDelta']))
    result.append(case)
    return result

def _tate_calculate_split(E, S, pAdics, T, case, result=[], **kwds):
    if case.has_key('roots'):
        if case['roots'] > 0:
            case['split'] = True
        else:
            case['split'] = False
        result.append(case)
        return result
    else:
        return get_number_of_roots_cases([S(1),
                                          S(E.a1()),
                                          S(-E.a2())],
                                         pAdics, T, 'roots', case,
                                         result=result, **kwds)
    
def _tate_finish(case, restrictions, result=[], variables=None, **kwds):
    tree = pAdicTree(variables=variables, root=case['T'])
    if len(restrictions) == 0:
        f = case['f']
        if f == 0:
            red_type = None
        elif f == 1:
            if case['split']:
                red_type = 1
            else:
                red_type = -1
        else:
            red_type = 0
        myresult = FreyCurveLocalData(case['E0'],
                                      case['T'].pAdics().prime_ideal(),
                                      f,
                                      case['vDelta'],
                                      KodairaSymbol(case['KS']),
                                      case['c'],
                                      red_type)
    else:
        myresult = []
        for r in restrictions:
            if r == 'conductor':
                myresult.append(case['f'])
            if r == 'reduction_type':
                f = case['f']
                if f == 0:
                    myresult.append(None)
                elif f == 1:
                    if case['split']:
                        myresult.append(1)
                    else:
                        myresult.append(-1)
                else:
                    myresult.append(0)
            if r == 'discriminant':
                myresult.append(case['vDelta'])
            if r == 'type':
                myresult.append(KodairaSymbol(case['KS']))
            if r == 'minimal_model':
                myresult.append(case['E0'])
    result.append((myresult, tree))
    return result
    
def _should_calculate_vDelta(case, restrictions):
    return not case.has_key('vDelta') and \
           ((len(restrictions) == 0) or
            ('discriminant' in restrictions) or
            (_should_calculate_m(case, restrictions)) or
            (_should_calculate_f(case, restrictions)) or
            (_should_calculate_n(case, restrictions)))

def _should_calculate_m(case, restrictions):
    return not case.has_key('m') and \
           ((len(restrictions) == 0) or
            (_should_calculate_f(case, restrictions)))

def _should_calculate_f(case, restrictions):
    return not case.has_key('f') and \
           ((len(restrictions) == 0) or
            ('conductor' in restrictions) or
            ('reduction_type' in restrictions))

def _should_calculate_n(case, restrictions):
    return case['KS'].count('n') > 0 and \
           ((len(restrictions) == 0) or
            ('type' in restrictions))

def _should_calculate_split(case, restrictions):
    return (case['f'] == 1 and
            not case.has_key('split') and
            (_should_calculate_c(case, restrictions) or
             'reduction_type' in restrictions))

def _should_calculate_c(case, restrictions):
    return not case.has_key('c') and \
           ((len(restrictions) == 0))
           
def _tate_cleanup(cases):
    result = []
    for case in cases:
        flag = True
        for case0 in result:
            if case0[0] == case[0]:
                # We can merge here to prevent unneccesary copying.
                # Since we would discard both old trees anyways this
                # is not a problem.
                case0[1].pAdic_tree()._root.merge(case[1]._root)
                flag = False
                break
        if flag:
            result.append((case[0], TreeCondition(case[1])))
    return ConditionalValue(result)

def _check_elliptic_curve(E):
   if not isinstance(E, EllipticCurve_generic):
            raise ValueError("%s is not an elliptic curve."%(E,))

def _init__coefficient_ring(coefficient_ring, E):
    if coefficient_ring is None:
        coefficient_ring = E.base_ring()
    for a in E.a_invariants():
        if a not in coefficient_ring:
            raise ValueError("%s is not part of %s."%(a, coefficient_ring))
    return coefficient_ring

def _init__pAdics(pAdics, ring, prime, coefficient_ring):
    if pAdics is None:
        if ring is None:
            ring = coefficient_ring.base_ring()
        if prime is None:
            raise ValueError("At least the argument pAdics or the argument prime should be given.")
        pAdics = pAdicBase(ring, prime)
    if not isinstance(pAdics, pAdicBase):
        raise ValueError("%s is not a pAdicBase object."%(pAdics,))
    return pAdics
 
def _init__polynomial_ring(coefficient_ring, pAdics):
    return coefficient_ring.change_ring(pAdics.number_field())
    
def _init__variables(polynomial_ring):
    return list(polynomial_ring.gens())
    
def _init__initial_values(initial_values, pAdics, variables):
    if initial_values is None:
        initial_values = pAdicTree(pAdics=pAdics, full=True,
                                   variables=variables)
    if not isinstance(initial_values, pAdicTree):
        raise ValueError("%s is not a p-adic tree."%(initial_value,))
    if not pAdics.is_extension_of(initial_values.pAdics()):
        raise ValueError("%s does not extend the p-adics of %s."%(pAdics,
                                                                  initial_values))
    return initial_values
    
def _init__cases(T, E):
    firstCase = dict(next_step=1, T=T.root(), E=E, E0=E)
    return [firstCase], []
    
def _init__str_list(str_list):
    if isinstance(str_list, str) or not hasattr(str_list, '__iter__'):
        str_list = [str_list]
    for s in str_list:
        if not isinstance(s, str):
            raise ValueError("%s is not a string."%(s,))
    return str_list
    
def performTatesAlgorithm(elliptic_curve, coefficient_ring=None, 
                          pAdics=None, base_ring=None, prime=None,
                          initial_values=None, verbose=False,
                          precision_cap=20, only_calculate=[]):
    r"""
    Performs Tate's Algorithm on an elliptic curve dependant on
    parameters.
    
    INPUT:
    
    - ``elliptic_curve`` -- An elliptic curve of which the
      coefficients of its Weierstrass equation are inside the
      polynomial ring coefficient_ring. Tate's algorithm will be
      performed on this elliptic curve.
      
    - ``coefficient_ring`` -- A polynomial ring
      (default: elliptic_curve.base_ring()) containing the
      coefficients of the elliptic curve. Note that this ring
      might be extended
      
    - ``pAdics`` -- A pAdicBase object
      (default: pAdicBase(base_ring, prime)) containing the
      prime and the associated number field over which we want
      to perform Tate's algorithm. Note that the arguments
      ``prime`` and ``base_ring`` will be ignored if this
      argument is given. This pAdicBase object must extend the
      p-adics of ``initial_values`` if it is given.
      
    - ``base_ring`` -- (default: coefficient_ring.base_ring()) A
      number field, or subring thereof, on which Tate's algorithm
      can be performed. If not specified it will be the coefficient
      ring of coefficient_ring which only works if coefficient_ring
      is actually a polynomial ring. In most cases this will be the
      ring of integers of a number field.
      
    - ``prime`` -- A prime (default: None) inside ``base_ring``
      for which Tate's algorithm must be performed. This may
      be given as an ideal in base_ring or a generator thereof.
      Note that this argument must be provided if the
      argument ``pAdics`` is not.
      
    - ``initial_values`` -- A pAdicTree (default: None) containing
      the possible values for the parameters in the given
      elliptic curve. If not given, this will be initialized by
      this function.
      
    - ``verbose`` -- A boolean value or an integer (default: False).
      When set to True or any value larger then zero will print
      comments to stdout about the computations being done whilst
      busy. If set to False or 0 will not print such comments. If
      set to any negative value will also prevent the printing of
      any warnings.
      If this method calls any method that accepts an argument
      verbose will pass this argument to it. If such a method
      fulfills a minor task within this method and the argument
      verbose was larger than 0, will instead pass 1 less than the
      given argument. This makes it so a higher value will print
      more details about the computation than a lower one.
      
    - ``precision_cap`` -- A non-negative integer (default: 20).
      This argument determines the highest precision that will be
      used for the approximation of the parameters. Note that
      setting this too low might result in inaccurate results, for
      which a warning will appear if this is the case. Setting
      this argument too high might result in endless and slow
      computations.
    
    - ``only_calculate`` -- (default: []) When set to a
      specific value or list of values, will ensure that the
      algorithm only calculates a certain quantity or quantities.
      Possible values are:
        - 'conductor' > Only calculate the conductor exponent
        - 'reduction_type' > Only calculate the reduction type, i.e.
          good, split multiplicative, non-split multiplicative or
          additive, returned as None, 1, -1 or 0 respectively.
        - 'discriminant' > Only calculate the valuation of the
                           discriminant
        - 'type' > Only calculate the Kodaira Symbol
        - 'minimal_model' > Only calculate the minimal model
                            for this elliptic curve
      Note that some calculations might require some other
      quantities to be determined as well, so the reduction in
      computation might depend on what choices are necessary.
    
    OUTPUT:
    
    The output is a list of possible outcomes of Tate's algorithm
    when performed on the elliptic curve ``elliptic_curve`` for
    the specified prime. When the argument only_calculate is not
    specified, each entry in this list is an object of the type
    ParametrizedLocalData containing the local information about
    this object, combined with the necessary values of the
    parameters stored in a pAdicTree. If however only_calculate
    was specified each entry will be a list wherein each entry
    corresponds to the requested value in only_calculate and
    the last entry is the pAdicTree of value for which these
    values are achieved.
    
    EXAMPLES:
    
    To be made
    
    """
    _check_elliptic_curve(elliptic_curve)
    coefficient_ring = _init__coefficient_ring(coefficient_ring, elliptic_curve)
    pAdics = _init__pAdics(pAdics, base_ring, prime, coefficient_ring)
    S = _init__polynomial_ring(coefficient_ring, pAdics)
    variables = _init__variables(S)
    T = _init__initial_values(initial_values, pAdics, variables)
    newCases, doneCases = _init__cases(T, elliptic_curve)
    only_calculate = _init__str_list(only_calculate)
    
    #The main loop performing the different steps.
    while len(newCases) > 0:
        cases = newCases
        newCases = []
        for case in cases:
            if case.has_key('next_step'):
                if case['next_step'] == 1:
                    if verbose > 0:
                        print "Performing Step 1"
                    _tate_step_1(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 2:         
                    if verbose > 0:
                        print "Performing Step 2 transformation"  
                    _tate_step_2_transformation(case['E'], S, pAdics, case['T'],
                                                case['E0'],
                                                variables=variables,
                                                result=newCases,
                                                verbose=(verbose-1 if verbose>0 else verbose),
                                                precision_cap=precision_cap)
                elif case['next_step'] == 2 + 1/2:
                    if verbose > 0:
                        print "Performing Step 2"
                    _tate_step_2(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 3:
                    if verbose > 0:
                        print "Performing Step 3"
                    _tate_step_3(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 4:
                    if verbose > 0:
                        print "Performing Step 4"
                    _tate_step_4(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 5:
                    if verbose > 0:
                        print "Performing Step 5"
                    _tate_step_5(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 6:
                    if verbose > 0:
                        print "Performing Step 6 transformation"
                    _tate_step_6_transformation(case['E'], S, pAdics, case['T'],
                                                case['E0'],
                                                variables=variables,
                                                result=newCases,
                                                verbose=(verbose-1 if verbose>0 else verbose),
                                                precision_cap=precision_cap)
                elif case['next_step'] == 6 + 1/2:
                    if verbose > 0:
                        print "Performing Step 6"
                    _tate_step_6(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 7:
                    if verbose > 0:
                        print "Performing Step 7"
                    _tate_step_7(case['E'], S, pAdics, case['T'], case['E0'],
                                 case, only_calculate, variables=variables,
                                 result=newCases, verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] in QQ and case['next_step'] < 8:
                    n , b = decodeQuotient( case['next_step'] - 7 )
                    if b == 0:
                        if verbose > 0:
                            print "Performing Step 7sub transformation"
                        _tate_step_7sub_transformation(case['E'], S, pAdics, case['T'],
                                                       case['E0'], n,
                                                       variables=variables,
                                                       result=newCases,
                                                       verbose=(verbose-1 if verbose>0 else verbose),
                                                       precision_cap=precision_cap)
                    else:
                        if verbose > 0:
                            print "Performing Step 7sub"
                        _tate_step_7sub(case['E'], S, pAdics, case['T'], case['E0'], n,
                                        variables=variables, result=newCases,
                                        verbose=(verbose-1 if verbose>0 else verbose),
                                        precision_cap=precision_cap)
                elif case['next_step'] == 8:
                    if verbose > 0:
                        print "Performing Step 8 transformation"
                    _tate_step_8_transformation(case['E'], S, pAdics, case['T'],
                                                case['E0'],
                                                variables=variables,
                                                result=newCases,
                                                verbose=(verbose-1 if verbose>0 else verbose),
                                                precision_cap=precision_cap)
                elif case['next_step'] == 8 + 1/2:
                    if verbose > 0:
                        print "Performing Step 8"
                    _tate_step_8(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 9:
                    if verbose > 0:
                        print "Performing Step 9 transformation"
                    _tate_step_9_transformation(case['E'], S, pAdics, case['T'],
                                                case['E0'],
                                                variables=variables,
                                                result=newCases,
                                                verbose=(verbose-1 if verbose>0 else verbose),
                                                precision_cap=precision_cap)
                elif case['next_step'] == 9 + 1/2:
                    if verbose > 0:
                        print "Performing Step 9"
                    _tate_step_9(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 10:
                    if verbose > 0:
                        print "Performing Step 10"
                    _tate_step_10(case['E'], S, pAdics, case['T'], case['E0'],
                                 variables=variables, result=newCases,
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap)
                elif case['next_step'] == 11:
                    if verbose > 0:
                        print "Performing Step 11"
                    _tate_step_11(case['E'], S, pAdics, case['T'], case['E0'],
                                  variables=variables, result=newCases,
                                  verbose=(verbose-1 if verbose>0 else verbose),
                                  precision_cap=precision_cap)
                else:
                    print "Unknown step number %s requested"%case['next_step']
            elif _should_calculate_vDelta(case, only_calculate):
                _tate_calculate_vDelta(case['E'], S, pAdics, case['T'], case,
                                       variables=variables, result=newCases,
                                       verbose=verbose,
                                       precision_cap=precision_cap)
            elif _should_calculate_m(case, only_calculate):
                _tate_calculate_m(case['E'], S, pAdics, case['T'], case,
                                  variables=variables, result=newCases,
                                  verbose=verbose, precision_cap=precision_cap)
            elif _should_calculate_f(case, only_calculate):
                _tate_calculate_f(case['E'], S, pAdics, case['T'], case,
                                  variables=variables, result=newCases,
                                  verbose=verbose, precision_cap=precision_cap)
            elif _should_calculate_n(case, only_calculate):
                _tate_calculate_n(case['E'], S, pAdics, case['T'], case,
                                  variables=variables, result=newCases,
                                  verbose=verbose, precision_cap=precision_cap)
            elif _should_calculate_split(case, only_calculate):
                _tate_calculate_split(case['E'], S, pAdics, case['T'], case,
                                      variables=variables, result=newCases,
                                      verbose=verbose, precision_cap=precision_cap)
            elif _should_calculate_c(case, only_calculate):
                _tate_calculate_c(case['E'], S, pAdics, case['T'], case,
                                  variables=variables, result=newCases,
                                  verbose=verbose, precision_cap=precision_cap)
            else:
                _tate_finish(case, only_calculate, result=doneCases,
                             verbose=verbose, variables=variables)
        
    return _tate_cleanup(doneCases)
    
def _calculate_vc4(E, S, pAdics, T, case, **kwds):
    return get_cases_invariant(S(E.c4()), pAdics, T, 'vc4', case, **kwds)
    
def _calculate_vc6(E, S, pAdics, T, case, **kwds):
    return get_cases_invariant(S(E.c6()), pAdics, T, 'vc6', case, **kwds)

def _papadopoulus_lookup_char_p(new_cases, S, variables=None, verbose=False,
                                precision_cap=20, restrictions=[]):
    result=[]
    while len(new_cases) > 0:
        cases = new_cases
        new_cases = []
        for case in cases:
            if not case.has_key('vDelta'):
                _tate_calculate_vDelta(case['E'], S, pAdics, case['T'], case,
                                       variables=variables, result=new_cases,
                                       verbose=verbose,
                                       precision_cap=precision_cap)
            elif case['vDelta'] == 0:
                case['KS'] = "I0"
                case['f'] = 0
                result.append(case)
            elif not case.has_key('vc4'):
                _calculate_vc4(case['E'], S, pAdics, case['T'], case,
                               variables=variables, result=new_cases,
                               verbose=verbose,
                               precision_cap=precision_cap)
            elif case['vc4'] == 0:
                n = case['vDelta']
                case['KS'] = "I"+str(n)
    
    return result

def perform_papadopoulus_lookup(elliptic_curve, initial_values=None,
    prime=None, coefficient_ring=None, base_ring=None, coPrimality=0,
    verbose=False, precision_cap=20, only_calculate=[]):
    r"""
    Calculates local data by doing a table lookup as provided in the
    article by Papadopoulus.
    
    INPUT:
    
    - ``elliptic_curve`` -- An elliptic curve of which the
      coefficients of its Weierstrass equation are inside the
      polynomial ring coefficient_ring. Tate's algorithm will be
      performed on this elliptic curve.
      
    - ``initial_values`` -- A pAdicTree (default: None) containing
      the possible values for the parameters in the given
      elliptic curve. If not given, this will be initialized by
      this function, in which case both the arguments prime and
      coefficient_ring must be defined.
      
    - ``prime`` -- A prime
      (default: initial_values.pAdics().prime_ideal()) inside
      ``base_ring`` for which Tate's algorithm must be performed.
      This may be given as an ideal in base_ring or a generator
      thereof.
      
    - ``coefficient_ring`` -- A polynomial ring over base_ring
      (default: initial_values.polynomial_ring()) of which the
      variables are precisely the parameters of elliptic_curve.
      
    - ``base_ring`` -- (default: coefficient_ring.base_ring()) A
      number field, or subring thereof, on which Tate's algorithm
      can be performed. If not specified it will be the coefficient
      ring of coefficient_ring which only works if coefficient_ring
      is actually a polynomial ring. In most cases this will be the
      ring of integers of a number field.
      
    - ``coPrimality`` -- A non-negative integer (default: 0).
      This option can be set to an integer n if it is known
      that the parameters are n-wise coprime, where n is the
      value of this parameter. The initial values for the
      parameters will then be set up such that no n parameters
      are simultaneously divisible by the given prime. If an
      initial case was set with the argument ``initial_values``
      then this argument is ignored. This argument can at most
      be the number of parameters.
      
    - ``verbose`` -- A boolean value (default: False). If set
      to True the program will print information about the steps
      it is performing and will ask its subfunctions to do the same.
      
    - ``precision_cap`` -- A non-negative integer (default: 20).
      This argument determines the highest precision that will be
      used for the approximation of the parameters. Note that
      setting this too low might result in inaccurate results, for
      which a warning will appear if this is the cases. Setting
      this argument too high might result in endless and slow
      computations.
    
    - ``only_calculate`` -- (default: []) When set to a
      specific value or list of values, will ensure that the
      algorithm only calculates a certain quantity or quantities.
      Possible values are:
        - 'conductor' > Only calculate the conductor exponent
        - 'discriminant' > Only calculate the valuation of the
                           discriminant
        - 'type' > Only calculate the Kodaira Symbol
        - 'minimal_model' > Only calculate the minimal model
                            for this elliptic curve
      Note that some calculations might require some other
      quantities to be determined as well, so the reduction in
      computation might depend on what choices are necessary.
    
    OUTPUT:
    
    The output is a list of possible outcomes of Tate's algorithm
    when performed on the elliptic curve ``elliptic_curve`` for
    the specified prime. When the argument only_calculate is not
    specified, each entry in this list is an object of the type
    ParametrizedLocalData containing the local information about
    this object, combined with the necessary values of the
    parameters stored in a pAdicTree. If however only_calculate
    was specified each entry will be a list wherein each entry
    corresponds to the requested value in only_calculate and
    the last entry is the pAdicTree of value for which these
    values are achieved.
    
    EXAMPLES:
    
    To be made
    
    """
    if initial_values is None:
        if coefficient_ring is None or prime is None:
            raise ValueError("If no initial values are given, both the arguments coefficient_ring and prime must be given.")
        if base_ring is None:
            base_ring = coefficient_ring.base_ring()
        variables = list(coefficient_ring.variable_names())
        pAdics = pAdicBase(base_ring, prime)
        initial_values = pAdicTree(root=pAdicNode(pAdics=pAdics, full=True,
                                                  width=len(variables)),
                                   variables=variables)
        if coPrimality not in ZZ or coPrimality < 0 or coPrimality > pAdics.width:
            raise ValueError("The argument coPrimality is not well-defined")
        if coPrimality > 0:
            initial_values.apply_coprimality_restriction(coPrimality=coPrimality)
    if isinstance(only_calculate, str):
        only_calculate = [only_calculate]
    
    firstCase = dict(next_step=1, T=initial_values.root(), E=elliptic_curve,
                     E0=elliptic_curve)
    cases = [firstCase]
    S = initial_values.polynomial_ring()
    variables = initial_values.variables()
    pAdics = initial_values.pAdics()
    char = pAdics.characteristic()
    if char == 2:
        l = pAdics.valuation(2)
        if l == 1:
            cases = _papadopoulus_lookup_char2_val1(cases, S, verbose=verbose,
                                                    variables=variables,
                                                    precision_cap=precision_cap,
                                                    restrictions=only_calculate)
        else:
            cases = _papadopoulus_lookup_char2(cases, S, l, verbose=verbose,
                                               precision_cap=precision_cap,
                                               variables=variables,
                                               restrictions=only_calculate)
    elif char == 3:
        l = pAdics.valuation(3)
        if l == 1:
            cases = _papadopoulus_lookup_char3_val1(cases, S, verbose=verbose,
                                                    precision_cap=precision_cap,
                                                    variables=variables,
                                                    restrictions=only_calculate)
        else:
            cases = _papadopoulus_lookup_char3(cases, S, l, verbose=verbose,
                                               precision_cap=precision_cap,
                                               variables=variables,
                                               restrictions=only_calculate)     
    else:
        cases = _papadopoulus_lookup_char_p(cases, S, verbose=verbose,
                                            precision_cap=precision_cap,
                                            variables=variables,
                                            restrictions=only_calculate)
    result = []
    for case in cases:
        _tate_finish(case, only_calculate, result=result,
                     verbose=verbose, variables=variables)
    return _tate_cleanup(result)
