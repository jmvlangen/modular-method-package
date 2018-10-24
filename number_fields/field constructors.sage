def field_with_root(K, a, names='sqrt_a', give_embedding=False):
    r"""
    Gets the field extension of K that contains the specified root.

    INPUT:
    
    - ``K`` -- A number field
    - ``a`` -- An element of the number field K
    - ``names`` -- A string (default: 'sqrt_a') or list thereof
      indicating the name(s) to use for the generators of the
      bigger field.
    - ``give_embedding`` -- A boolean (default: False) that
      indicates whether embeddings of K into this field and
      a square root of should be returned.
    
    OUPUT:

    A number field that contains K and a square root of a. If it is
    an extension of K and not K itself its generators will have the
    names specified by names. If give_embeddings was set to True, will
    return a tuple consisting of
     - the number field as mentioned before
     - an embedding of K into that number field
     - a square root of a inside that number field

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
         2/5*a^3 - 3/5*a^2 + 2*a - 2/5)
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
            K_to_L = K.hom([a.minpoly().change_ring(L).roots()[0][0] for a in K.gens()], L)
            return L, K_to_L, sqrt(K_to_L(a))
        else:
            return L

def fixed_field(H):
    r"""
    Returns the fixed field of a subset of a galois group

    INPUT:

    - ``H`` -- An iterable object containing elements of
      a galois group. len(H) should be at least 1

    OUTPUT:

    A number field K consisting of all those elements that
    are mapped to themselves by elements of H.

    EXAMPLES:

    A simple example::

        sage: K = CyclotomicField(12)
        sage: G = K.galois_group()
        sage: H = [G.gens()[0]]
        sage: fixed_field(H)
        Number Field in zeta120 with defining polynomial x^2 - 2*x + 4

    If H only contains the trivial element, the entire
    field is returned::

        sage: K = CyclotomicField(12)
        sage: G = K.galois_group()
        sage: H = [G.identity()]
        sage: fixed_field(H)
        Cyclotomic Field of order 12 and degree 4

    H empty does not work::

        sage: fixed_field([])
        Traceback (most recent call last)
        ...
        IndexError: list index out of range

    If H generates or is the entire galois group we get
    the rational field::

        sage: K = CyclotomicField(24)
        sage: G = K.galois_group()
        sage: fixed_field(G)
        Rational Field
        sage: fixed_field(G.gens())
        Rational Field
    """
    G = H[0].parent()
    if H == G:
        return QQ
    if hasattr(H, 'fixed_field'):
        result = H.fixed_field()
        if isinstance(result, tuple):
            return result[0]
        else:
            return result
    return fixed_field(G.subgroup(H))

@cached_function
def composite_field(K1, K2, give_maps=False):
    r"""
    Returns the composite field of K1 and K2

    INPUT:

    - ``K1`` -- A number field
    - ``K2`` -- A number field
    - ``give_maps`` -- A boolean (default=False) indicating whether
      the embeddings should be returned.

    OUTPUT:

    A number field K that is the composite field of K1 and K2. If
    give_maps was set to True, will instead return a tuple consisting of
     - The field K
     - An embedding of K1 into K
     - An embedding of K2 into K

    EXAMPLES:

    Combining two quadratic fields::

        sage: K1 = QuadraticField(2)
        sage: K2 = QuadraticField(3)
        sage: K = composite_field(K1, K2); K
        Number Field in a0 with defining polynomial x^4 - 10*x^2 + 1
        sage: K(2).is_square() and K(3).is_square()
        True

    Also works if one of the fields contains the other::

        sage: K1 = QuadraticField(2)
        sage: K2 = CyclotomicField(8)
        sage: K = composite_field(K1, K2); K
        Cyclotomic Field of order 8 and degree 4
        sage: K2(2).is_square()
        True

    Can use the optional give_maps to obtain the embeddings::

        sage: K1 = QuadraticField(2)
        sage: K2 = QuadraticField(3)
        sage: composite_field(K1, K2, give_maps=True)
        (Number Field in a0 with defining polynomial x^4 - 10*x^2 + 1, Ring morphism:
           From: Number Field in a with defining polynomial x^2 - 2
           To:   Number Field in a0 with defining polynomial x^4 - 10*x^2 + 1
           Defn: a |--> -1/2*a0^3 + 9/2*a0, Ring morphism:
           From: Number Field in a with defining polynomial x^2 - 3
           To:   Number Field in a0 with defining polynomial x^4 - 10*x^2 + 1
           Defn: a |--> -1/2*a0^3 + 11/2*a0)
    """
    if not is_field(K1):
        if K1.is_subring(QQ):
            K1 = QQ
        else:
            K1 = K1.field_of_fractions()
    if not is_field(K2):
        if K2.is_subring(QQ):
            K2 = QQ
        else:
            K2 = K2.field_of_fractions()
    from_K2 = None
    if K1 != QQ and K2 != QQ and K1.defining_polynomial().parent() != K2.defining_polynomial().parent():
        R = K1.defining_polynomial().parent()
        f2 = R(K2.defining_polynomial())
        K2orig = K2
        K2 = NumberField(f2, names=K2.variable_names())
        if give_maps:
            from_K2 = K2orig.hom(K2.gens(), K2)
    if K1.is_subring(K2):
        if give_maps:
            if from_K2 is None:
                return K2, K1.embeddings(K2)[0], K2.embeddings(K2)[0]
            else:
                return K2, K1.embeddings(K2)[0], K2.embeddings(K2)[0] * from_K2
        else:
            return K2
    elif K2.is_subring(K1):
        if give_maps:
            if from_K2 is None:
                return K1, K1.embeddings(K1)[0], K2.embeddings(K1)[0]
            else:
                return K1, K1.embeddings(K1)[0], K2.embeddings(K1)[0] * from_K2
        else:
            return K1
    else:
        if give_maps:
            if from_K2 is None:
                return K1.composite_fields(K2, both_maps=give_maps)[0][0:3]
            else:
                result = list(K1.composite_fields(K2, both_maps=give_maps)[0][0:3])
                result[2] = result[2] * from_K2
                return tuple(result)
        else:
            return K1.composite_fields(K2)[0]
