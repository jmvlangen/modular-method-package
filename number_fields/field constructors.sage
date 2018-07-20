def field_with_root(K, a, names='sqrt_a'):
    r"""
    Gets the field extension of K that contains the specified root.

    INPUT:
    
    - ``K`` -- A number field
    - ``a`` -- An element of the number field K
    - ``names`` -- A string (default: 'sqrt_a') or list thereof
                   indicating the name(s) to use for the generators
                   of the bigger field.
    
    OUPUT:

    A number field that contains K and a square root of a. If it is
    an extension of K and not K itself its generators will have the
    names specified by names.
    """
    a = K(a)
    if a.is_square():
        return K
    else:
        R.<x> = K[]
        return K.extension(x^2 - a, names=names)

def fixed_field(H):
    r"""
    Returns the fixed field of a subset of a galois group

    INPUT:

    - ``H`` -- An iterable object containing elements of
               a galois group. len(H) should be at least 1

    OUTPUT:

    A number field K consisting of all those elements that
    are mapped to themselves by elements of H.
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
    """
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
