def galois_field_extend(sigma, K, embedding=None):
    r"""
    Finds the extension of a galois homomorphism to a bigger field.

    INPUT:

    - ``sigma`` -- A galois homomorphism on some field K0
    - ``K`` -- A galois field extension of K0
    - ``embedding`` -- An embedding of K0 into K. By default
                       will be the map sending K0.gen() to
                       the first found root of its minimal
                       polynomial in K

    OUTPUT:

    An element tau of the galois group of K, such that tau
    restricted to K0 is sigma.
    """
    Gsmall = sigma.parent()
    Ksmall = Gsmall.number_field()
    Gbig = K.galois_group()
    Kbig = Gbig.number_field()
    if Gsmall == Gbig and Ksmall == Kbig:
        return sigma
    if embedding is None:
        embedding = Ksmall.hom([a.minpoly().change_ring(Kbig).roots()[0][0] for a in Ksmall.gens()], Kbig)
    for tau in Gbig:
        if tau(embedding(Ksmall.gen())) == embedding(sigma(Ksmall.gen())):
            return tau
    raise Exception("No corresponding galois action found for %s."%sigma)

def galois_field_restrict(sigma, K, embedding=None):
    r"""
    Finds the restriction of a galois homomorphism to a smaller field.

    INPUT:

    - ``sigma`` -- A galois homomorphism on some field K0
    - ``K`` -- A galois subfield of K0
    - ``embedding`` -- An embedding of K into K0. By default
                       will be the map sending K.gen() to the
                       first found root of its minimal polynomial
                       in K0.

    OUPUT:
    
    An element tau of the galois group of K such that sigma and tau
    acts the same on elements of K.
    """
    Gsmall = K.galois_group()
    Ksmall = Gsmall.number_field()
    Gbig = sigma.parent()
    Kbig = Gbig.number_field()
    if Gbig == Gsmall and Ksmall == Kbig:
        return sigma
    if embedding is None:
        embedding = Ksmall.hom([a.minpoly().change_ring(Kbig).roots()[0][0] for a in Ksmall.gens()], Kbig)
    for tau in Gsmall:
        if sigma(embedding(Ksmall.gen())) == embedding(tau(Ksmall.gen())):
            return tau
    raise Exception("No corresponding galois action found for %s."%tau)

def galois_field_change(sigma, K):
    r"""
    Changes a Galois homomorphism to one on a specified field
    
    Given a galois homomorphism sigma on a field K0, will first
    find an extension of sigma to K0 K and then restrict that
    extension to K.

    INPUT:
    
    - ``sigma`` -- An element of the galois group of some field K0
    - ``K`` -- A galois field extension of the same base field as K0

    OUTPUT:

    An element tau of the galois group of K, such that sigma and tau
    have a common extension to the compositum of K0 and K, i.e. there
    is some galois homomorphism mu on K0 K that acts the same as
    sigma on K0 and the same as tau on K. Note that if K0 and K have
    no common subfield, sigma and tau might not have anything in common.
    """
    L = sigma.parent().number_field()
    M, L_to_M, K_to_M = composite_field(L, K, give_maps=True)
    return galois_field_restrict(galois_field_extend(sigma, M, embedding=L_to_M), K, embedding=K_to_M)

def cyclotomic_galois_isomorphism(s, N=None):
    r"""
    Realizes the isomorphism between the galois group of a cyclotomic field and $\Z/N\Z^*$

    This function works both ways, converting integers into galois elements
    and vice versa.

    INPUT:

    - ``s`` -- An integer with gcd(s,N) == 1 or an element
               of the galois group of a number field.
    - ``N`` -- A strictly positive integer indicating the
               order of the cyclotomic field. Should be specified
               if s is an integer. Otherwise defaults to the
               size of the torsion of the corresponding number
               field.
    
    OUTPUT:

    If s was an integer, will return an element \sigma of the galois
        group of $\Q(\zeta_N)$ such that $ \sigma(\zeta_N) = \zeta_N^s $
    If s was an element of the galois group of a number field K, will
        return an integer n such that s(zeta) == zeta^n, where zeta
        is a generator of the torsion of K. If N was specified will
        do the same but with the field $\Q(\zeta_N)$ instead of
        K and instead of s a galois homomorphism on $\Q(\zeta_N)$
        that has a common extension to the composite field of K and
        $\Q(\zeta_N)$
    """
    if s in ZZ:
        if N is None:
            raise Exception("Need to specify N if s is an integer")
        L.<zeta> = CyclotomicField(N)
        G = L.galois_group()
        for sigma in G:
            if sigma(zeta) == zeta^s:
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
    i = 0
    zeta = L.zeta(N)
    zeta_s = s(zeta)
    while zeta_s != 1:
        zeta_s = zeta_s/zeta
        i += 1
    return i
