r"""Code to work with galois homomorphisms

This code provides methods to change the definition field of a galois
homomorphism. Besides that it also contains code to convert galois
homomorphisms of $\Q(\zeta_N)$ into the corresponding elements of
$(\Z/N\Z)^*$ and back. There is also a method to see how a galois
group acts on the units of the field.

EXAMPLES:

We can extend and restrict galois homomorphisms in towers of number
fields, but the construction only works well in one direction::

    sage: from modular_method.number_fields.galois_group import galois_field_extend, galois_field_restrict
    sage: K1 = QuadraticField(3)
    sage: K2 = CyclotomicField(12)
    sage: G1.<s> = K1.galois_group()
    sage: G2.<t,u> = K2.galois_group()
    sage: galois_field_extend(s, K2) == t
    True
    sage: galois_field_restrict(t, K1) == s
    True
    sage: galois_field_restrict(u, K1) == s
    True

The method :func:`galois_field_change` finds the best way to convert
between two fields::

    sage: from modular_method.number_fields.galois_group import galois_field_change
    sage: K1 = QuadraticField(3)
    sage: K2 = CyclotomicField(12)
    sage: K3 = CyclotomicField(20)
    sage: G1.<s> = K1.galois_group()
    sage: G2.<t,u> = K2.galois_group()
    sage: G3.<v,w> = K3.galois_group()
    sage: galois_field_change(s, K2) == t
    True
    sage: galois_field_change(v, K2) == u^4
    True
    sage: galois_field_change(u, K3) == v^4
    True

For cyclotomic fields we can jump back and forth::

    sage: from modular_method.number_fields.galois_group import cyclotomic_galois_isomorphism
    sage: N = 120
    sage: n = 19
    sage: s = cyclotomic_galois_isomorphism(n, N=N); s
    (1,15)(2,16)(3,13)(4,14)(5,11)(6,12)(7,9)(8,10)(17,31)(18,32)(19,29)(20,30)(21,27)(22,28)(23,25)(24,26)
    sage: cyclotomic_galois_isomorphism(s, N=N)
    19

We can inspect how the galois group acts on the units of a number
field::

    sage: from modular_method.number_fields.galois_group import galois_on_units
    sage: K = QuadraticField(3)
    sage: G.<s> = K.galois_group()
    sage: U = K.unit_group()
    sage: d = galois_on_units(K)
    sage: exp = [1, -2]
    sage: U(s(product(U.gens()[i]^exp[i] for i in range(2))))
    u0*u1^2
    sage: exps = vector(exp)*d[s]
    sage: product(U.gens()[i]^exps[i] for i in range(2))
    u0*u1^2

AUTHORS:

- Joey van Langen (2019-02-15): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_function

from sage.all import ZZ, QQ, Integer

from sage.matrix.constructor import matrix

from sage.rings.number_field.number_field import CyclotomicField
from modular_method.number_fields.field_constructors import intersection_field, composite_field

@cached_function(key=lambda s, K, e: (str(s), s.parent().number_field(), K, e))
def galois_field_extend(sigma, K, embedding=None):
    r"""Find an extension of a galois homomorphism to a bigger field.

    INPUT:

    - ``sigma`` -- A galois homomorphism on some field $K_0$

    - ``K`` -- A galois field extension of $K_0$

    - ``embedding`` -- An embedding of $K_0$ into `K`. If the Galois group
      of `K` is abelian will default to the first embedding found by
      the method :meth:`embeddings`. Otherwise will raise an error if
      not provided.

    OUTPUT:

    An element $\tau$ of the galois group of `K`, such that $\tau$
    restricted to $K_0$ along the given embedding is `sigma`.

    EXAMPLES::

    Simple example::

        sage: from modular_method.number_fields.galois_group import galois_field_extend, galois_field_restrict
        sage: K = QuadraticField(2)
        sage: L = CyclotomicField(8)
        sage: sigma = K.galois_group().gen()
        sage: tau = galois_field_extend(sigma, L); tau
        (1,3)(2,4)
        sage: tau.parent()
        Galois group of Cyclotomic Field of order 8 and degree 4
        sage: galois_field_restrict(tau, K) == sigma
        True

    Note that the extension is not unique::

        sage: from modular_method.number_fields.galois_group import galois_field_extend, galois_field_restrict
        sage: K = QuadraticField(2)
        sage: L = CyclotomicField(8)
        sage: tau = L.galois_group()[3]
        sage: sigma = galois_field_restrict(tau, K)
        sage: galois_field_extend(sigma, L) == tau
        False

    If the smaller field is not abelian, the extension actually
    depends on the chosen embedding. Therefore an error is raised if
    none is given::

        sage: from modular_method.number_fields.galois_group import galois_field_extend, galois_field_restrict
        sage: R.<x> = QQ[]
        sage: K0.<a0> = NumberField(x^3 - 2)
        sage: K.<a> = K0.galois_closure()
        sage: L0.<b0> = NumberField(x^6 - 2)
        sage: L.<b> = L0.galois_closure()
        sage: sigma = K.galois_group().gens()[0]
        sage: galois_field_extend(sigma, L)
        Traceback (most recent call last):
        ...
        ValueError: The extension of Galois elements is not well-defined if no embedding is given.
        sage: tau = [galois_field_extend(sigma, L, embedding=phi) for phi in K.embeddings(L)]
        sage: all(tau[i] == tau[0] for i in range(len(tau)))
        False
        sage: all(galois_field_restrict(tau[i], K, embedding=K.embeddings(L)[0]) == sigma
        ....:     for i in range(len(tau)))
        False

    """
    Gsmall = sigma.parent()
    Ksmall = Gsmall.number_field()
    Gbig = K.galois_group()
    Kbig = Gbig.number_field()
    if Gsmall == Gbig and Ksmall == Kbig:
        return sigma
    if embedding is None:
        if not Gsmall.is_abelian():
            raise ValueError("The extension of Galois elements is not " +
                             "well-defined if no embedding is given.")
        embedding = Ksmall.embeddings(Kbig)[Integer(0) ]
    for tau in Gbig:
        if tau(embedding(Ksmall.gen())) == embedding(sigma(Ksmall.gen())):
            return tau
    raise Exception("No corresponding galois action found for %s."%sigma)

@cached_function(key=lambda s, K, e: (str(s), s.parent().number_field(), K, e))
def galois_field_restrict(sigma, K, embedding=None):
    r"""Find the restriction of a galois homomorphism to a smaller field.

    INPUT:

    - ``sigma`` -- A galois homomorphism on some field $K_0$

    - ``K`` -- A galois subfield of $K_0$

    - ``embedding`` -- An embedding of K into K0. If the Galois group
      of `K` is abelian will default to the first embedding found by
      the method :meth:`embeddings`. Otherwise will raise an error if
      not provided.

    OUPUT:
    
    An element $\tau$ of the galois group of `K` such that `sigma` and
    $\tau$ act the same on elements of `K`.

    EXAMPLES::

        sage: from modular_method.number_fields.galois_group import galois_field_restrict
        sage: K = QuadraticField(3)
        sage: L = CyclotomicField(24)
        sage: tau = L.galois_group().gens()[0]
        sage: sigma = galois_field_restrict(tau, K); sigma
        (1,2)
        sage: sigma.parent()
        Galois group of Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?

    If the smaller field is not abelian, the restriction is not unique::

        sage: from modular_method.number_fields.galois_group import galois_field_restrict
        sage: R.<x> = QQ[]
        sage: K0.<a0> = NumberField(x^3 - 2)
        sage: K.<a> = K0.galois_closure()
        sage: L0.<b0> = NumberField(x^6 - 2)
        sage: L.<b> = L0.galois_closure()
        sage: sigma = L.galois_group().gens()[0]
        sage: galois_field_restrict(sigma, K)
        Traceback (most recent call last):
        ...
        ValueError: The restriction of Galois elements is not unique if no embedding is given.
        sage: tau = [galois_field_restrict(sigma, K, embedding=phi) for phi in K.embeddings(L)]
        sage: all(tau[i] == tau[0] for i in range(len(tau)))
        False

    """
    Gsmall = K.galois_group()
    Ksmall = Gsmall.number_field()
    Gbig = sigma.parent()
    Kbig = Gbig.number_field()
    if Gbig == Gsmall and Ksmall == Kbig:
        return sigma
    if embedding is None:
        if not Gsmall.is_abelian():
            raise ValueError("The restriction of Galois elements is not " +
                             "unique if no embedding is given.")
        embedding = Ksmall.embeddings(Kbig)[Integer(0) ]
    for tau in Gsmall:
        if sigma(embedding(Ksmall.gen())) == embedding(tau(Ksmall.gen())):
            return tau
    raise Exception("No corresponding galois action found for %s."%sigma)

@cached_function(key=lambda s, K, L: (str(s), s.parent().number_field(), K, L))
def galois_field_change(sigma, K, L=None):
    r"""Change a Galois homomorphism to one on a specified field
    
    Given a galois homomorphism $\sigma$ and a galois number field
    $K$, will find a galois homomorphism $\tau$ of $K$ such that
    $\sigma$ and $\tau$ have a common extension to a galois
    homomorphism of the algebraic closure of $\Q$.

    .. NOTE::

    The galois homomorphism returned by this method only agrees with
    the given galois homomorphism on the intersection of their
    respective number fields.

    If the intersection of these two number fields is not abelian, the
    restriction to this intersection depends on the embeddings. Hence
    if no specific embeddings are given, this function raises an
    error.

    INPUT:
    
    - ``sigma`` -- An element of the galois group of some number field
      $K_0$.

    - ``K`` -- A galois number field.

    - ``L`` -- A tuple consisting of a number field containing both
      $K_0$ and `K`, an embedding from $K_0$ into that field, and an
      emedding from `K` into that field, in that order. If set to
      `None` (default) it will be initialized as a composite field
      using the function :func:`composite_field` with the
      corresponding embeddings.

    OUTPUT:

    An element $\tau$ of the galois group of `K`, such that `sigma`
    and $\tau$ have a common extension to the compositum of $K_0$ and
    `K`, i.e. there is some galois homomorphism $\mu$ on $K_0 K$ that
    acts the same as `sigma` on $K_0$ and the same as $\tau$ on
    `K`.

    EXAMPLES:

    Can convert from a very complicated field to a very simple one::

        sage: from modular_method.number_fields.galois_group import galois_field_change
        sage: K = QuadraticField(10)
        sage: L = CyclotomicField(120)
        sage: tau = L.galois_group()[5]; tau
        (1,6,3,8)(2,7,4,5)(9,14,11,16)(10,15,12,13)(17,22,19,24)(18,23,20,21)(25,30,27,32)(26,31,28,29)
        sage: sigma = galois_field_change(tau, K); sigma
        (1,2)
        sage: sigma.parent()
        Galois group of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380?

    Can also be used on two totally unrelated fields, but the result
    can be totally unexpected::

        sage: from modular_method.number_fields.galois_group import galois_field_change
        sage: K = QuadraticField(2)
        sage: L = QuadraticField(3)
        sage: sigma = K.galois_group().gen(); sigma
        (1,2)
        sage: tau = galois_field_change(sigma, L); tau
        ()
        sage: tau.parent()
        Galois group of Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?

    Note that if the intersection of the two fields is not abelian,
    the change in Galois homomorphism depends on chosen embeddings::

        sage: from modular_method.number_fields.galois_group import galois_field_change
        sage: R.<x> = QQ[]
        sage: K1.<a1> = NumberField(x^6 - 2)
        sage: K1.<a1> = K1.galois_closure()
        sage: K2.<a2> = NumberField(x^6 - 3*x^4 - 2*x^3 + 3*x^2 + 6*x + 1)
        sage: K2.<a2> = K2.galois_closure()
        sage: G.<sigma, tau> = K1.galois_group()
        sage: galois_field_change(sigma, K2)
        Traceback (most recent call last):
        ...
        ValueError: The intersection field must be abelian to make restrictions of Galois homomorphism unique.
        sage: from modular_method.number_fields.field_constructors import composite_field
        sage: L = composite_field(K1, K2)
        sage: galois_field_change(sigma, K2, L=(L, K1.embeddings(L)[0], K2.embeddings(L)[0]))
        (1,9)(2,10)(3,11)(4,12)(5,6)(7,8)
        sage: galois_field_change(sigma, K2, L=(L, K1.embeddings(L)[0], K2.embeddings(L)[1]))
        (1,10)(2,6)(3,8)(4,11)(5,12)(7,9)

    .. SEEALSO::

        :func:`galois_field_restrict`,
        :func:`galois_field_extend`,
        :func:`composite_field`,
        :func:`intersection_field`,

    TESTS:

    No unnecessary errors are raised when the fields are equal, even
    when this field is not abelian::

        sage: from modular_method.number_fields.galois_group import galois_field_change
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^6 + 3)
        sage: G = K.galois_group()
        sage: s = G[2]
        sage: s == galois_field_change(s, K)
        True

    """
    K1 = sigma.parent().number_field()
    if K1 == K:
        return sigma
    K0, K0_to_K1, K0_to_K = intersection_field(K1, K, give_maps=True, L=L)
    if K0 == QQ:
        return K.galois_group()[Integer(0) ]
    if not K0.is_abelian() and L is None:
        raise ValueError("The intersection field must be abelian to make " +
                         "restrictions of Galois homomorphism unique.")
    return galois_field_extend(galois_field_restrict(sigma, K0,
                                                     embedding=K0_to_K1),
                               K, embedding=K0_to_K)

@cached_function(key=lambda s, N: ((Integer(0) , s, N) if s in ZZ else
                                   (str(s), s.parent().number_field(), N)))
def cyclotomic_galois_isomorphism(s, N=None):
    r"""Realize the isomorphism between the galois group of a cyclotomic
    field and $\Z/N\Z^*$

    There is a natural isomorphism between the group $(\Z/N\Z)^*) and
    the galois group of $\Q(\zeta_N)$ where the integer $n$
    corresponds to the galois homomorphism $\zeta_N \mapsto
    \zeta_N^n$. This function allows one to convert the element of
    one of these groups into the other.

    INPUT:

    - ``s`` -- An integer with gcd(s,N) == 1 or an element of the
      galois group of a number field.

    - ``N`` -- A strictly positive integer indicating the order of the
      cyclotomic field. Should be specified if s is an
      integer. Otherwise defaults to the size of the torsion of the
      corresponding number field.
    
    OUTPUT:

    If `s` was an integer, will return an element $\sigma$ of the galois
    group of $\Q(\zeta_N)$ such that $ \sigma(\zeta_N) = \zeta_N^s $

    If `s` was an element of the galois group of a number field $K$,
    will return an integer `n` such that $s(\zeta) == \zeta^n$, where
    $\zeta$ is a generator of the roots of unity in $K$. If `N` was
    specified will do the same but will rather replace $K$ by the
    field $\Q(\zeta_N)$, replace `s` by a galois homomorphism of that
    $\Q(\zeta_N)$ that has a common extension with `s`, and replace
    $\zeta$ by $\zeta_N$.

    EXAMPLES:

    Converting from an integer to a galois homomorphism::

        sage: from modular_method.number_fields.galois_group import cyclotomic_galois_isomorphism
        sage: sigma = cyclotomic_galois_isomorphism(3, N=8); sigma
        (1,4)(2,3)
        sage: sigma.parent()
        Galois group of Cyclotomic Field of order 8 and degree 4

    Converting from a galois homomorphism to the corresponding
    integer::

        sage: from modular_method.number_fields.galois_group import cyclotomic_galois_isomorphism
        sage: L = CyclotomicField(30)
        sage: sigma = L.galois_group()[2]; sigma
        (1,3)(2,4)(5,7)(6,8)
        sage: cyclotomic_galois_isomorphism(sigma)
        19

    If the modulus is given can convert from a galois homomorphism
    over any field::

        sage: from modular_method.number_fields.galois_group import cyclotomic_galois_isomorphism
        sage: K = QuadraticField(5)
        sage: sigma = K.galois_group().gen()
        sage: cyclotomic_galois_isomorphism(sigma, N=5)
        2

    """
    if s in ZZ:
        if N is None:
            raise Exception("Need to specify N if s is an integer")
        L = CyclotomicField(N, names=('zeta',)); (zeta,) = L._first_ngens(1)
        G = L.galois_group()
        for sigma in G:
            if sigma(zeta) == zeta**s:
                return sigma
        raise Exception("No galois element corresponds to %s."%s)
    G = s.parent()
    L = G.number_field()
    if N is None:
        N = L.zeta_order()
    else:
        s = galois_field_change(s, CyclotomicField(N))
        G = s.parent()
        L = G.number_field()
    i = Integer(0) 
    zeta = L.zeta(N)
    zeta_s = s(zeta)
    while zeta_s != Integer(1) :
        zeta_s = zeta_s/zeta
        i += Integer(1) 
    return i

@cached_function
def galois_on_units(K):
    r"""Compute how the galois group of a field acts on its units.

    INPUT:

    - ``K`` -- A galois number field.
    
    OUTPUT:
    
    A dictionary indexed by the elements of the galois group of `K`,
    such that the value with key $\sigma$ is a matrix with integer
    coefficients. If $u_0$, ..., $u_n$ are the generators of the unit
    group of `K` as given by :meth:`unit_group` and $a_0$, ..., $a_n$
    are the entries of the i-th row of this matrix, then $u_i$ is
    mapped to $u_0^{a_0} * ... * u_n^{a_n} by s.

    EXAMPLE::

        sage: from modular_method.number_fields.galois_group import galois_on_units
        sage: K = QuadraticField(3)
        sage: G.<s> = K.galois_group()
        sage: U = K.unit_group()
        sage: d = galois_on_units(K)
        sage: exp = [1, -2]
        sage: U(s(product(U.gens()[i]^exp[i] for i in range(2))))
        u0*u1^2
        sage: exps = vector(exp)*d[s]
        sage: product(U.gens()[i]^exps[i] for i in range(2))
        u0*u1^2

    """
    G = K.galois_group()
    U = K.unit_group()
    return {s : matrix([U(s(u)).list() for u in U.gens()]).transpose() for s in G}

