r"""
Implements twisting methods for elliptic curves

This file contains special methods that can twist a curve
and find all possible twists that change the conductor.

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

def _is_twist_good_char(E1, E2):
    """
    Determines whether the curve E1 is a twist of E2.

    INPUT:

    - ``E1`` -- An elliptic curve over a field of
      characteristic not 2 or 3.
    - ``E2`` -- An ellitpic curve defined over the
      same field as E1.
    
    OUTPUT:
    1 if the two curves are isomorphic
    -1 if the curves are twists of one another,
    but not isomorphic over their base field.
    0 if the curves are not twists of one another.
    """
    if E1.c4()^3 * E2.c6()^2 != E2.c4()^3 * E1.c6()^2:
        return 0
    if E2.c4() == 0:
        c = E1.c6() / E2.c6()
        try:
            a = c.nth_root(3)
            if a.is_square():
                return 1
            else:
                return -1
        except ValueError:
            return 0
    elif E2.c6() == 0:
        b = E1.c4() / E2.c4()
        try:
            a = b.nth_root(2)
            if a.is_square():
                return 1
            else:
                return -1
        except ValueError:
            return 0
    else:
        b = E1.c4() / E2.c4()
        c = E1.c6() / E2.c6()
        a = c / b
        if a.is_square():
            return 1
        else:
            return -1

def _is_twist_char_2(E1, E2):
    """
    Determines whether the curve E1 is a twist of E2.

    INPUT:

    - ``E1`` -- An elliptic curve over a field of
      characteristic 2.
    - ``E2`` -- An ellitpic curve defined over the
      same field as E1.
    
    OUTPUT:
    1 if the two curves are isomorphic
    -1 if the curves are twists of one another,
    but not isomorphic over their base field.
    0 if the curves are not twists of one another.
    """
    
    if (E1.a1() != 0 or
        E1.a3() != 0 or
        E2.a1() != 0 or
        E2.a3() != 0):
        return 0
    b1 = E1.a2()^2 + E1.a4()
    b2 = E2.a2()^2 + E2.a4()
    c1 = E1.a2()*E1.a4() + E1.a6()
    c2 = E2.a2()*E2.a4() + E2.a6()
    if b1^3 * c2^2 != b2^3 * c1^2:
        return 0
    if b2 == 0:
        c = c1 / c2
        try:
            a = c.nth_root(3)
            if a.is_square():
                return 1
            else:
                return -1
        except ValueError:
            return 0
    elif c2 == 0:
        b = b1 / b2
        try:
            a = b.nth_root(2)
            if a.is_square():
                return 1
            else:
                return -1
        except ValueError:
            return 0
    else:
        b = b1 / b2
        c = c1 / c2
        a = c / b
        if a.is_square():
            return 1
        else:
            return -1    

def _is_twist_char_3(E1, E2):
    """
    Determines whether the curve E1 is a twist of E2.

    INPUT:

    - ``E1`` -- An elliptic curve over a field of
      characteristic 3.
    - ``E2`` -- An ellitpic curve defined over the
      same field as E1.
    
    OUTPUT:
    1 if the two curves are isomorphic
    -1 if the curves are twists of one another,
    but not isomorphic over their base field.
    0 if the curves are not twists of one another.
    """
    if (E1.b2()^2 * E2.b4() != E2.b2()^2 * E1.b4() or
        E1.b2()^3 * E2.b6() != E2.b2()^2 * E1.b6() or
        E1.b4()^3 * E2.b4()^2 != E2.b4()^3 * E1.b6()^2):
        return 0
    if E2.b2() == 0:
        if E2.b4() == 0:
            c = E1.b6() / E2.b6()
            try:
                a = c.nth_root(3)
                if a.is_square():
                    return 1
                else:
                    return -1
            except ValueError:
                return 0
        elif E2.b6() == 0:
            b = E1.b4() / E2.b4()
            try:
                a = b.nth_root(2)
                if a.is_square():
                    return 1
                else:
                    return -1
            except ValueError:
                return 0
        else:
            b = E1.b4() / E2.b4()
            c = E1.b6() / E2.b6()
            a = c / b
            if a.is_square():
                return 1
            else:
                return -1
    else:
        a = E1.b2() / E2.b2()
        if a.is_square():
            return 1
        else:
            return -1

def is_twist(E1, E2):
    """
    Determines whether the curve E1 is a twist of E2.

    INPUT:

    - ``E1`` -- An elliptic curve over a field.
    - ``E2`` -- An ellitpic curve defined over the
      same field as E1.
    
    OUTPUT:
    1 if the two curves are isomorphic
    -1 if the curves are twists of one another,
    but not isomorphic over their base field.
    0 if the curves are not twists of one another.
    """
    p = E1.base_ring().characteristic()
    if p == 2:
        return _is_twist_char_2(E1, E2)
    elif p == 3:
        return _is_twist_char_3(E1, E2)
    else:
        return _is_twist_good_char(E1, E2)

def twist_elliptic_curve(E, d):
    """
    Returns the twists of an elliptic curve by 'd'.

    Twisting an elliptic curve with Weierstrass equation
    ..MATH::
    
        y^2 = x^3 + a_2 x^2 + a_4 x + a_6

    means changing it into the curve given by
    ..MATH::
    
        y^2 = x^3 + d a_2 x^2 + d^2 a_4 x + d^3 a_6

    which is isomorphic to the first curve over $R(\\sqrt{d})$,
    where $R$ is the ring over which the curve was defined.

    INPUT:

    - ``E`` -- An elliptic curve given by a Weierstrass equation
               with $a_1 = 0$ and $a_3 = 0$
    - ``d`` -- A ring element that can be multiplied with the
               coefficients of E. In general this will be an element
               of the base ring of E, but any ring for which coercion
               with those elements is defined should work.

    OUTPUT:
    
    An elliptic curve that is the twist of the given curve E by the
    parameter d.

    EXAMPLES:
    
    A simple examples ::

        sage: E = EllipticCurve([1,2]); E
        Elliptic Curve defined by y^2 = x^3 + x + 2 over Rational Field
        sage: twist_elliptic_curve(E,-1)
        Elliptic Curve defined by y^2 = x^3 + x - 2 over Rational Field

    An example over a more complicated field ::

        sage: K = CyclotomicField(5)
        sage: K.<zeta> = CyclotomicField(5)
        sage: E = EllipticCurve([zeta,3]); E
        Elliptic Curve defined by y^2 = x^3 + zeta*x + 3 over Cyclotomic Field of order 5 and degree 4
        sage: twist_elliptic_curve(E, zeta)
        Elliptic Curve defined by y^2 = x^3 + zeta^3*x + 3*zeta^3 over Cyclotomic Field of order 5 and degree 4

    We can also twist using parameters ::

        sage: R.<a> = QQ[]
        sage: E = EllipticCurve([a^2 + 3, a-1]); E
        Elliptic Curve defined by y^2 = x^3 + (a^2+3)*x + (a-1) over Univariate Polynomial Ring in a over Rational Field
        sage: twist_elliptic_curve(E, a+1)
        Elliptic Curve defined by y^2 = x^3 + (a^4+2*a^3+4*a^2+6*a+3)*x + (a^4+2*a^3-2*a-1) over Univariate Polynomial Ring in a over Rational Field
    """
    if E.a1() != 0 or E.a3() != 0:
        raise ValueError("Can only twist if a1 and a3 are zero.")
    a2 = d   * E.a2()
    a4 = d^2 * E.a4()
    a6 = d^3 * E.a6()
    return EllipticCurve([0,a2,0,a4,a6])

def compute_possible_twists(K, P):
    r"""
    Computes twists that could change the conductor of an elliptic curve.

    Given an elliptic curve over a number field, the conductor exponent
    at a given prime might change if the curve is twisted. This function
    computes a list of elements of the number field, such that all these
    changes in the conductor at a given prime can be achieved by twisting
    by one of the elements in this list.

    INPUT:
    
    - ``K`` -- A number field, possibly $\\Q$, over which the elliptic
               curve would be defined.
    - ``P`` -- A prime ideal of K.

    OUTPUT:
    A list of elements of K. Every change in the conductor exponent
    at P of an elliptic curve defined over K, corresponds to at least
    one twist of the curve by an element in this list. The trivial
    twist (1) is always part of this list.

    EXAMPLES:

    Well known results for $\\Q$ ::

        sage: compute_possible_twists(QQ, 2)
        [1, 3, 2, 6]
        sage: compute_possible_twists(QQ, 5)
        [1, 5]

    Works over bigger number fields ::

        sage: K = CyclotomicField(7)
        sage: compute_possible_twists(K, K.prime_above(2))
        [1,
        -zeta7^5 + zeta7^4 - zeta7^3,
        -zeta7^5 + zeta7^4 + zeta7^3,
        -zeta7^5 - zeta7^4 + zeta7^3,
        zeta7^5 - zeta7^4 + zeta7^3,
        zeta7^5 + zeta7^4 - zeta7^3,
        -zeta7^5 - zeta7^4 - zeta7^3,
        zeta7^5 + zeta7^4 + zeta7^3,
        zeta7^3 + zeta7 + 1,
        2*zeta7^5 + 2*zeta7^4 + zeta7^3 + 2*zeta7^2 + zeta7 + 3,
        2*zeta7^4 + zeta7^3 - zeta7 + 1,
        -2*zeta7^5 + zeta7^3 - zeta7 - 1,
        -2*zeta7^5 - 2*zeta7^4 - zeta7^3 - 2*zeta7^2 - zeta7 - 3,
        2*zeta7^5 - zeta7^3 + zeta7 + 1,
        zeta7^3 + 2*zeta7^2 + zeta7 + 1,
        -zeta7^3 - 2*zeta7^2 - zeta7 - 1]
        sage: compute_possible_twists(K, K.prime_above(3))
        [1, 3]

    """
    base = pAdicBase(K, P)
    pi = base.uniformizer()
    if base.characteristic() != 2:
        return [K(1), pi]
    else:
        result =  _compute_possible_unit_twists(base)
        return result + [c * pi for c in result]

def _compute_possible_unit_twists(base):
    r"""
    Returns the different twists of an elliptic curve by units.

    INPUT:
    - ``base`` -- A pAdicBase object that contains the field and prime

    OUTPUT:
    Similar to the output of :func: compute_possible_twists,
    but only computes those that are units modulo the prime.
    """
    K = base.number_field()
    T_unram = _compute_unramified_unit_tree(base)
    T_ram = T_unram.complement()
    T_ram.root().children.get((K(0),)).remove()
    k = T_unram.root().minimum_full_level()
    unram = [node.quotient_tuple()[0] for node in T_unram.nodes_at_level(k)]
    ram = [node.quotient_tuple()[0] for node in T_ram.nodes_at_level(k)]
    result = [K(1)]
    while len(ram) > 0:
        val = ram[0]
        result.append(val.lift())
        for val2 in unram:
            if val * val2 in ram:
                ram.remove(val * val2)
            else:
                "Warning: possibly a mistake in `_compute_possible_unit_twists`"
    return result

def _tree_of_unit_squares_mod(base, mod):
    r"""
    Compute the squares modulo some power of a prime in a local ring.

    INPUT:
    - ``base`` -- a pAdicBase object that describes the local ring.
    - ``mod`` -- the power of the prime modulo which we want to
                 consider squares.

    OUTPUT:
    A pAdicTree that contains all numbers in the local ring,
    that are squares modulo the prime to the power mod.
    """
    K = base.number_field()
    k = base.valuation(2)
    n = min(mod - k, floor((mod + 1)/2))
    T = pAdicTree(variables='x', pAdics=base)
    T.root().children.get((K(0),)).remove()
    T2 = pAdicTree(variables='x', pAdics=base, full=False)
    for node in T.nodes_at_level(n):
        val = node.representative()[0]^2
        node = T2.root()
        for coeff in base.power_series(val, mod):
            coeffs = (coeff,)
            if not node.children.contains(coeffs):
                node.children.add(pAdicNode(parent=node,
                                            coefficients=coeffs))
            node = node.children.get(coeffs)
        if not node.is_full():
            node.children = pAdicNodeCollection_inverted(node)
    return T2
        
def _compute_unramified_unit_tree(base):
    r"""
    Computes the tree of units of a local ring that give unramified
    quadratic extensions.

    INPUT:
    - ``base`` -- A pAdicBase object that describes the local ring

    OUTPUT:
    A pAdicTree that contains all units of the local ring, such that
    the extension by the square root of these remains unramified.
    """
    k = base.valuation(2)
    T = _tree_of_unit_squares_mod(base, 2*k)
    for i in reversed(range(0,2*k,2)):
        for node in T.nodes_at_level(i):
            for coeffs in base.representatives():
                if not node.children.contains(coeffs) \
                and not (i == 0 and coeffs[0] == 0):
                    node.children.add(pAdicNode(parent=node,
                                                coefficients=coeffs,
                                                full=True))
    return T
