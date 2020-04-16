r"""Ways of constructing fields.

Some special methods to construct fields that have uniform output.
These methods include a function to construct the field extension
created by adding a root, a function to construct the composite field
of two fields and a function 

EXAMPLES:

We can build field extensions by adding roots::

    sage: K1.<a> = field_with_root(QQ, 2); K1
    Number Field in a with defining polynomial x^2 - 2
    sage: K2.<b> = field_with_root(K1, 3); K2
    Number Field in b with defining polynomial x^4 - 10*x^2 + 1
    sage: K3.<c> = field_with_root(K2, 6); K3
    Number Field in b with defining polynomial x^4 - 10*x^2 + 1

We can find the fixed field of some galois elements::

    sage: L = CyclotomicField(24)
    sage: G = L.galois_group()
    sage: fixed_field([G[1], G[4]])
    Number Field in zeta240 with defining polynomial x^2 - 2
    sage: fixed_field([G[5]])
    Number Field in zeta240 with defining polynomial x^4 - 4*x^2 + 1

We can find the composite field of two number fields, even if one of
them is $\Q$::

    sage: K1.<a> = QuadraticField(-7)
    sage: K2.<b> = CyclotomicField(24)
    sage: K3.<c> = QQ[sqrt(2), sqrt(3)]
    sage: composite_field(K1, K2)
    Number Field in ab with defining polynomial x^16 + 56*x^14 + 1370*x^12 + 19236*x^10 + 169739*x^8 + 960036*x^6 + 3382958*x^4 + 6637820*x^2 + 5536609
    sage: composite_field(K1, K3)
    Number Field in asqrt2 with defining polynomial x^8 + 8*x^6 + 256*x^4 + 3648*x^2 + 14400
    sage: composite_field(K2, K3)
    Cyclotomic Field of order 24 and degree 8
    sage: composite_field(QQ, K3)
    Number Field in xsqrt2 with defining polynomial x^4 - 10*x^2 + 1
    sage: composite_field(QQ, K2)
    Cyclotomic Field of order 24 and degree 8

AUTHORS:

- Joey van Langen (2019-02-15): initial version

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

def field_with_root(K, a, names='sqrt_a', give_embedding=False):
    r"""Get a field extension of K that contains the specified root.

    INPUT:
    
    - ``K`` -- A number field, which may be the rationals.

    - ``a`` -- An element of the number field K.

    - ``names`` -- A string (default: 'sqrt_a') or list thereof
      indicating the name(s) to use for the generators of the bigger
      field.

    - ``give_embedding`` -- A boolean (default: False) that indicates
      whether embeddings of K into this field and a square root of `a`
      should be returned.
    
    OUPUT:

    A number field that contains `K` and a square root of `a`. If it
    is an extension of `K` and not `K` itself its generators will have
    the names specified by names. If give_embeddings was set to True,
    will return a tuple consisting of: the number field as mentioned
    before; an embedding of K into that number field; and a square
    root of `a` inside that number field

    EXAMPLES:

    Simple examples over $\Q$::
    
        sage: field_with_root(QQ, 3)
        Number Field in sqrt_a with defining polynomial x^2 - 3
        sage: field_with_root(QQ, -2)
        Number Field in sqrt_a with defining polynomial x^2 + 2
    
    Working over a bigger field also works::

        sage: K = CyclotomicField(5); K
        Cyclotomic Field of order 5 and degree 4
        sage: field_with_root(K, -2, names='a')
        Number Field in a with defining polynomial x^8 - 2*x^7 + 11*x^6 - 16*x^5 + 39*x^4 - 28*x^3 + 19*x^2 + 6*x + 11

    The root might already be contained in the field::

        sage: K = CyclotomicField(5); K
        Cyclotomic Field of order 5 and degree 4
        sage: field_with_root(K, 5)
        Cyclotomic Field of order 5 and degree 4
        sage: K(5).is_square()
        True

    A map can also be generated::

        sage: K = CyclotomicField(3); K
        Cyclotomic Field of order 3 and degree 2
        sage: field_with_root(K, -2, names='a', give_embedding=True)
        (Number Field in a with defining polynomial x^4 - 2*x^3 + 7*x^2 - 6*x + 3,
         Ring morphism:
           From: Cyclotomic Field of order 3 and degree 2
           To:   Number Field in a with defining polynomial x^4 - 2*x^3 + 7*x^2 - 6*x + 3
           Defn: zeta3 |--> 2/5*a^3 - 3/5*a^2 + 2*a - 7/5,
         2/5*a^3 - 3/5*a^2 + 3*a - 7/5)

    """
    a = K(a)
    if a.is_square():
        if give_embedding:
            return K, K.hom(K), sqrt(a)
        else:
            return K
    else:
        R.<x> = K[]
        L = K.extension(x^2 - a, names=names).absolute_field(names=names)
        if give_embedding:
            K_to_L = K.embeddings(L)[0]
            return L, K_to_L, sqrt(K_to_L(a))
        else:
            return L

def fixed_field(H, map=False):
    r"""Return the fixed field of a subset of a galois group

    INPUT:

    - ``H`` -- An iterable object containing elements of a galois
      group. len(H) should be at least 1

    - ``map`` -- A boolean value (default: False) indicating whether a
      map to the big field should be returned.

    OUTPUT:

    A number field K consisting of all those elements that are mapped
    to themselves by elements of H. It the argument `map` is set to
    True, will return a tuple of which the first part is the mentioned
    field K and the second part is an embedding of K into the galois
    number field corresponding to elements of H.

    EXAMPLES:

    A simple example::

        sage: K = CyclotomicField(12)
        sage: G = K.galois_group()
        sage: H = [G.gens()[0]]
        sage: fixed_field(H)
        Number Field in zeta120 with defining polynomial x^2 - 2*x + 4

    If H only contains the trivial element, the entire field is
    returned::

        sage: K = CyclotomicField(12)
        sage: G = K.galois_group()
        sage: H = [G.identity()]
        sage: fixed_field(H)
        Number Field in zeta120 with defining polynomial x^4 - x^2 + 1

    H empty does not work::

        sage: fixed_field([])
        Traceback (most recent call last)
        ...
        IndexError: list index out of range

    If H generates or is the entire galois group we get the rational
    field::

        sage: K = CyclotomicField(24)
        sage: G = K.galois_group()
        sage: fixed_field(G)
        Rational Field
        sage: fixed_field(G.gens())
        Rational Field

    Maps can be returned in each of the above cases::

        sage: K = CyclotomicField(12)
        sage: G = K.galois_group()
        sage: H = [G.identity()]
        sage: fixed_field([G.gens()[0]], map=True)
        (Number Field in zeta120 with defining polynomial x^2 - 2*x + 4, Ring morphism:
           From: Number Field in zeta120 with defining polynomial x^2 - 2*x + 4
           To:   Cyclotomic Field of order 12 and degree 4
           Defn: zeta120 |--> 2*zeta12^2)
        sage: fixed_field([G.identity()], map=True)
        (Cyclotomic Field of order 12 and degree 4,
         Identity endomorphism of Cyclotomic Field of order 12 and degree 4)
        sage: fixed_field(G, map=True)
        (Rational Field, Coercion map:
           From: Rational Field
           To:   Cyclotomic Field of order 12 and degree 4)

    """
    if hasattr(H, '_ambient'):
        G = H._ambient
    else:
        G = H[0].parent()
    if H == G:
        if map:
            return QQ, QQ.hom(G.number_field())
        else:
            return QQ
    if hasattr(H, 'fixed_field'):
        result = H.fixed_field()
        if isinstance(result, tuple):
            if map:
                return result
            else:
                return result[0]
        else:
            if map:
                return result, result.hom(G.number_field())
            else:
                return result
    return fixed_field(G.subgroup(H), map=map)

def write_as_extension(phi, give_map=False, names=None):
    r"""Give a field embedding as a field extension

    INPUT:

    - ``phi`` -- An embedding of fields

    - ``give_map`` -- A boolean value (default: False). Returns the
      map from the codomain of `phi` to the returned value as a second
      argument if set to True.

    - ``names`` -- A string with the name of the variable for the
      given field.

    OUTPUT:

    A number field $K$ that is an extension of the number field that
    is the domain of `phi` and isomorphic to the codomain of `phi`. If
    the argument `give_map` is set to True, will also return an
    isomorphism from the codomain of `phi` to the field $K$, such that
    this map combined with `phi` is precisely the coercion of the
    domain of `phi` into $K$.

    EXAMPLES::

        sage: write_as_extension(QuadraticField(5).embeddings(CyclotomicField(5))[0])
        Number Field in zeta5 with defining polynomial x^2 + (1/2*a + 1/2)*x + 1 over its base field
    
    One can name the variables as usual::

        sage: K.<a> = QuadraticField(2)
        sage: L.<b> = CyclotomicField(8)
        sage: M.<c> = write_as_extension(K.embeddings(L)[0])
        sage: M.absolute_polynomial() == L.absolute_polynomial()
        True
        sage: M.base_field() == K
        True

    One can also request the corresponding isomorphism. Note that its
    combination with the embeddings is simply the coercion map::

        sage: K.<a> = QuadraticField(3)
        sage: L.<b> = CyclotomicField(12)
        sage: M, phi = write_as_extension(K.embeddings(L)[0], give_map=True)
        sage: phi
        Ring morphism:
          From: Cyclotomic Field of order 12 and degree 4
          To:   Number Field in b with defining polynomial x^2 + a*x + 1 over its base field
          Defn: b |--> b
        sage: phi * K.embeddings(L)[0]
        Ring morphism:
          From: Number Field in a with defining polynomial x^2 - 3
          To:   Number Field in b with defining polynomial x^2 + a*x + 1 over its base field
          Defn: a |--> a

    TESTS::

    Still works if the codomain is the rationals::

        sage: write_as_extension(QQ.hom(QQ))
        Number Field in x1 with defining polynomial x - 1

    Also works if the codomain is an extension of fields::

        sage: K.<a> = QQ[sqrt(2), sqrt(3)]
        sage: write_as_extension(QQ.hom(K), names='b')
        Number Field in b with defining polynomial x^4 - 10*x^2 + 1

    Resolves conflicts in naming::

        sage: K.<a> = QuadraticField(2)
        sage: L.<a> = CyclotomicField(8)
        sage: write_as_extension(K.embeddings(L)[0])
        Number Field in a1 with defining polynomial x^2 + a*x + 1 over its base field
        sage: K.<a1> = QuadraticField(3)
        sage: L.<a1> = CyclotomicField(12)
        sage: write_as_extension(K.embeddings(L)[0])
        Number Field in a2 with defining polynomial x^2 + a1*x + 1 over its base field

    Works if both fields are the same::

        sage: K = QuadraticField(3)
        sage: write_as_extension(K.hom(K))
        Number Field in a1 with defining polynomial x - a over its base field

    Resolves conflicts in naming if a name was provided::

        sage: K.<a> = QuadraticField(2)
        sage: L.<a> = CyclotomicField(8)
        sage: write_as_extension(K.embeddings(L)[0], names='a')
        Number Field in a1 with defining polynomial x^2 + a*x + 1 over its base field

    """
    phi = _write_as_im_gen_map(phi)
    K = phi.domain()
    L = phi.codomain()
    if not is_field(K) and not is_field(L):
        raise ValueError("Input must be an embedding of fields.")
    if not L.is_absolute():
        L = L.absolute_field(names=L.variable_name())
        phi = _concat_maps(phi, L.structure()[1])
        concat_map = True
    else:
        concat_map = False
    if L == QQ:
        R.<x> = QQ[]
        f = x - 1
    else:
        f = L.absolute_polynomial()
    for g in f.change_ring(K).factor():
        if g[0].change_ring(phi)(L.gen()) == 0:
            if names is None:
                names = L.variable_name()
            if K.variable_name() == names:
                for n in range(len(names)):
                    if names[n:].isdigit():
                        names = names[:n] + str(Integer(names[n:]) + 1)
                        break
                else:
                    names = names + "1"
            M = K.extension(g[0], names=names)
            if give_map:
                phi = L.hom([M.gen()])
                if concat_map:
                    phi = _concat_maps(L.structure()[1], phi)
                return M, phi
            else:
                return M

@cached_function
def composite_field(K1, K2, give_maps=False, names=None):
    r"""Return the composite field of K1 and K2

    INPUT:

    - ``K1`` -- A number field, which may be the rationals. It may
      also be the extension of some number field $K_0$ of which `K2`
      is also an extension. This extension may also be given as an
      embedding of $K_0$ into `K1`.

    - ``K2`` -- A number field, which may be the rationals. It may
      also be the extension of some number field $K_0$ of which `K1`
      is also an extension. This extension may also be given as an
      embedding of $K_0$ into `K2`.

    - ``give_maps`` -- A boolean (default=False) indicating whether
      the embeddings into the composite field should be returned.

    - ``names`` -- A string, tuple thereof or None (default:
      `None`). If not `None` this will be used as the variable names
      in the composite field.

    OUTPUT:

    A number field $K$ that is the composite field of `K1` and
    `K2`. If `K1` and `K2` are given as extensions of a common field
    $K_0$ the field `K` will also be an extension of $K_0$ respecting
    the common structure. If `give_maps` was set to True, will instead
    return a tuple consisting of: the field K; an embedding of K1 into
    K; and an embedding of K2 into K. If one of the given fields `K1`
    and `K2` is absolute and isomorphic to the field $K$, then $K$
    will always be that field. In case both fields are absolute and
    isomorphic, then $K$ will be the field `K2`.

    EXAMPLES:

    Combining two quadratic fields::

        sage: K1 = QuadraticField(2)
        sage: K2 = QuadraticField(3)
        sage: K = composite_field(K1, K2); K
        Number Field in a1 with defining polynomial x^4 - 10*x^2 + 1
        sage: K(2).is_square() and K(3).is_square()
        True

    If one of the fields is contained in the other, returns the bigger
    field::

        sage: K1 = QuadraticField(2)
        sage: K2 = CyclotomicField(8)
        sage: K = composite_field(K1, K2)
        sage: K2 == K
        True
        sage: K2(2).is_square()
        True

    Note that if both fields are isomorphic, will always return the
    second field::

        sage: K1 = QuadraticField(-1)
        sage: K2 = QuadraticField(-4)
        sage: K1.is_isomorphic(K2)
        True
        sage: K = composite_field(K1, K2)
        sage: K == K2
        True
        sage: K == K1
        False

    Can use the optional give_maps to obtain the embeddings::

        sage: K1 = QuadraticField(2)
        sage: K2 = QuadraticField(3)
        sage: composite_field(K1, K2, give_maps=True)
        (Number Field in a1 with defining polynomial x^4 - 10*x^2 + 1, Ring morphism:
           From: Number Field in a with defining polynomial x^2 - 2
           To:   Number Field in a1 with defining polynomial x^4 - 10*x^2 + 1
           Defn: a |--> -1/2*a1^3 + 9/2*a1, Ring morphism:
           From: Number Field in a with defining polynomial x^2 - 3
           To:   Number Field in a1 with defining polynomial x^4 - 10*x^2 + 1
           Defn: a |--> -1/2*a1^3 + 11/2*a1)

    """
    to_K1 = None; to_K2 = None
    K1or = K1; K2or = K2
    if hasattr(K1, "codomain"):
        K1or = K1.codomain()
        K1 = write_as_extension(K1, give_map=give_maps)
        if give_maps:
            K1, to_K1 = K1
    if hasattr(K2, "codomain"):
        K2or = K2.codomain()
        K2 = write_as_extension(K2, give_map=give_maps)
        if give_maps:
            K2, to_K2 = K2
    if give_maps and to_K1 is None:
        to_K1 = K1.hom(K1)
    if give_maps and to_K2 is None:
        to_K2 = K2.hom(K2)
    if K2.absolute_degree() < K1.absolute_degree():
        if give_maps:
            K, K2_to_K, K1_to_K = composite_field(K2, K1,
                                                  give_maps=True)
            return (K,
                    _concat_maps(to_K1, K1_to_K),
                    _concat_maps(to_K2, K2_to_K))
        else:
            return composite_field(K2, K1)
    K01 = K1
    K02 = K2
    while K01 != K02:
        if K01.absolute_degree() < K02.absolute_degree():
            K02 = K02.base_ring()
        else:
            K01 = K01.base_ring()
    K1 = write_as_extension(K01.hom(K1), give_map=give_maps)
    K2 = write_as_extension(K02.hom(K2), give_map=give_maps)
    if give_maps:
        to_K1 = _concat_maps(to_K1, K1[1])
        K1 = K1[0]
        to_K2 = _concat_maps(to_K2, K2[1])
        K2 = K2[0]
    f = K2.relative_polynomial()
    g = f.change_ring(K1).factor()[0][0]
    if g.degree() * K1.absolute_degree() == K2.absolute_degree():
        if give_maps:
            if K01 == QQ:
                a0 = 1
                i1 = K01.hom(K1or)
                i2 = K01.hom(K2or)
            else:
                a0 = K01.absolute_generator()
                i1 = [i for i in K01.embeddings(K1or) if to_K1(i(a0)) == a0][0]
                i2 = [i for i in K01.embeddings(K2or) if to_K2(i(a0)) == a0][0]
            iota = [i for i in K1or.embeddings(K2or) if i(i1(a0)) == i2(a0)][0]
            if K2or.is_absolute():
                return K2or, iota, K2or.hom(K2or)
            else:
                K = K2or.absolute_field(names=K2or.variable_name())
                iota2 = _write_as_im_gen_map(K.structure()[1])
                iota = _concat_maps(iota, iota2)
                return K, iota, iota2
        else:
            if K2or.is_absolute():
                return K2or
            else:
                return K2or.absolute_field(names=K2or.variable_name())
    name1 = K1.variable_name(); name2 = K2.variable_name()
    if names is None:
        names = name1 + (name2 if name2 != name1 else "")
    if names == name1:
        for n in range(len(names)):
            if names[n:].isdigit():
                names = names[:n] + str(Integer(names[n:]) + 1)
                break
        else:
            names = names + "1"
    K = K1.extension(g, names=names)
    if give_maps:
        K1_to_K = _write_as_im_gen_map(K1.hom(K))
        K2_to_K = K2.hom([K.gen()])
        if not K.is_absolute():
            K = K.absolute_field(names=K.variable_name())
            K1_to_K = _concat_maps(K1_to_K, K.structure()[1])
            K2_to_K = _concat_maps(K2_to_K, K.structure()[1])
        return (K,
                _concat_maps(to_K1, K1_to_K),
                _concat_maps(to_K2, K2_to_K))
    else:
        if K.is_absolute():
            return K
        else:
            return K.absolute_field(names=K.variable_name())

@cached_function
def intersection_field(K1, K2, L=None, give_maps=False, names=None):
    r"""Give the intersection of two fields

    Given two number fields gives the intersection of these fields in
    a common field that contains both. This intersection is itself a
    field.

    INPUT:

    - ``K1`` -- A number field, which may be the rationals.
    
    - ``K2`` -- A number field, which may be the rationals.

    - ``L`` -- A tuple consisting of a number field containing both
      `K1` and `K2`, an embedding from `K1` to that field, and an
      embedding from `K2` to that field, in that order. If set to None
      (default) it will be initialized as a composite field using
      :func:`composite_field` with the corresponding embeddings.

    - ``give_maps`` -- A boolean (default: `False`) indicating whether
      the embeddings of the intersection field to `K1` and `K2` should
      be returned.

    - ``names`` -- A string, list thereof or None (default:
      `None`). If not `None` this will be used as the variable names
      in the intersection field.

    OUTPUT:

    A number field `K` isomorphic to the intersection of `K1` and `K2`
    in the number field `L`. If `give_maps` was set to `True`, will
    instead return a tuple consisting of: the field `K`; an embedding
    of `K` into `K1`; and an embedding of `K` into `K2`. These
    embeddings respect the embeddings given in `L`.

    EXAMPLES:

    A simple example::

        sage: R.<x> = QQ[]
        sage: K1.<a1> = CyclotomicField(8)
        sage: K2.<a2> = NumberField(x^4 - 2)
        sage: K0.<a0> = intersection_field(K1, K2); K0
        Number Field in a0 with defining polynomial x^2 - 2

    Using the method give_maps, one can obtain the corresponding
    embeddings::

        sage: R.<x> = QQ[]
        sage: K1.<a1> = CyclotomicField(9)
        sage: K2.<a2> = NumberField(x^4 + 3)
        sage: intersection_field(K1, K2, give_maps=True, names=a0)
        (Number Field in a0 with defining polynomial x^2 + 3, Ring morphism:
           From: Number Field in a0 with defining polynomial x^2 - x + 1
           To:   Cyclotomic Field of order 9 and degree 6
           Defn: a0 |--> -a1^3, Ring morphism:
           From: Number Field in a0 with defining polynomial x^2 - x + 1
           To:   Number Field in a2 with defining polynomial x^4 + 3
           Defn: a0 |--> -1/2*a2^2 + 1/2)

    Note that the intersection field might depend on the choice of
    common field::

        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^4 - x + 1)
        sage: L.<a> = K.galois_closure()
        sage: K1.<a1> = L.subfields(degree=12)[0][0]
        sage: K2.<a2> = [k for k, _, _ in L.subfields(degree=12) if not K1.is_isomorphic(k)][0]
        sage: intersection_field(K1, K2, L=(L, K1.embeddings(L)[0], K2.embeddings(L)[0]))
        Number Field in a0a6 with defining polynomial x^6 - x^4 - x^3 - x^2 + 1
        sage: intersection_field(K1, K2, L=(L, K1.embeddings(L)[1], K2.embeddings(L)[0]))
        Number Field in a0a6 with defining polynomial x^3 - 4*x - 1    

    .. SEEALSO::

        :func:`composite_field`

    """
    if K1.absolute_degree() > K2.absolute_degree():
        if give_maps:
            K, K_to_K2, K_to_K1 = intersection_field(K2, K1, L=L,
                                                     give_maps=True,
                                                     names=names)
            return K, K_to_K1, K_to_K2
        else:
            return intersection_field(K2, K1, L=L, names=names)
    if K1 == QQ:
        if give_maps:
            return QQ, QQ.hom(K1), QQ.hom(K2)
        else:
            return QQ
    if L is None:
        L, K1_to_L, K2_to_L = composite_field(K1, K2, give_maps=True,
                                              names=names)
    else:
        L, K1_to_L, K2_to_L = L
    if names is None:
        name1 = K1.variable_name(); name2 = K2.variable_name()
        names = name1 + (name2 if name2 != name1 else "")
    a1 = K1.absolute_generator()
    a2 = K2.absolute_generator()
    n1 = K1.absolute_degree()
    n2 = K2.absolute_degree()
    B1 = vector([a1^n for n in range(n1)])
    B2 = vector([a2^n for n in range(n2)])
    M1 = matrix([K1_to_L(b).list() for b in B1])
    M2 = matrix([K2_to_L(b).list() for b in B2])
    B = block_matrix([[M1], [-M2]]).kernel()
    n = B.dimension()
    if n == 1:
        if give_maps:
            return QQ, QQ.hom(K1), QQ.hom(K2)
        else:
            return QQ
    b1 = vector([B1 * b[0:n1] for b in B.basis()])
    t = 0
    while True:
        t += 1
        for ls in Partitions(t+n, length=n):
            for v in Permutations(ls):
                v = vector(list(v)) - vector(len(v)*[1])
                alpha1 = b1 * v
                f = alpha1.absolute_minpoly()
                if f.degree() == n:
                    K = NumberField(f, names=names)
                    Kopt, from_opt, to_opt = K.optimized_representation(name=names)
                    Kopt = Kopt.change_names(K.variable_name())
                    from_opt = from_opt * Kopt.structure()[0]
                    to_opt = Kopt.structure()[1] * to_opt
                    if give_maps:
                        b2 = vector([B2 * b[n1:] for b in B.basis()])
                        alpha2 = b2 * v
                        return (Kopt,
                                _concat_maps(from_opt, K.hom([alpha1])),
                                _concat_maps(from_opt, K.hom([alpha2])))
                    else:
                        return Kopt

def _write_as_im_gen_map(phi, dom=None):
    r"""Rewrites the map phi such that it is given by images of generators"""
    if dom is None:
        dom = phi.domain()
    g = dom.gen()
    if dom.is_absolute():
        return dom.hom([phi(g)])
    else:
        base_hom = _write_as_im_gen_map(phi, dom=dom.base_ring())
        F = dom.Hom(phi.codomain())
        return F(phi(g), base_map=base_hom)
                    
def _concat_maps(phi, psi):
    r"""Concatenates phi and psi, phi first"""
    return _write_as_im_gen_map(psi * phi)
