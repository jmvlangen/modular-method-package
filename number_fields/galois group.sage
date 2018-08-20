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

    EXAMPLES::

    Simple example::

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

        sage: K = QuadraticField(2)
        sage: L = CyclotomicField(8)
        sage: tau = L.galois_group()[3]
        sage: sigma = galois_field_restrict(tau, K)
        sage: galois_field_extend(sigma, L) == tau
        False
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

    EXAMPLE::

        sage: K = QuadraticField(3)
        sage: L = CyclotomicField(24)
        sage: tau = L.galois_group().gens()[0]
        sage: sigma = galois_field_restrict(tau, K); sigma
        (1,2)
        sage: sigma.parent()
        Galois group of Number Field in a with defining polynomial x^2 - 3
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

    EXAMPLES:

    Can convert from a very complicated field to a very simple one::

        sage: K = QuadraticField(10)
        sage: L = CyclotomicField(120)
        sage: tau = L.galois_group()[5]; tau
        (1,6,3,8)(2,7,4,5)(9,14,11,16)(10,15,12,13)(17,22,19,24)(18,23,20,21)(25,30,27,32)(26,31,28,29)
        sage: sigma = galois_field_change(tau, K); sigma
        (1,2)
        sage: sigma.parent()
        Galois group of Number Field in a with defining polynomial x^2 - 10

    Can also be used on two totally unrelated fields, but the result
    can be totally unexpected::

        sage: K = QuadraticField(2)
        sage: L = QuadraticField(3)
        sage: sigma = K.galois_group().gen(); sigma
        (1,2)
        sage: tau = galois_field_change(sigma, L); tau
        ()
        sage: tau.parent()
        Galois group of Number Field in a with defining polynomial x^2 - 3
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

    EXAMPLES:

    Converting from an integer to a galois homomorphism::

        sage: sigma = cyclotomic_galois_isomorphism(3, N=8); sigma
        (1,4)(2,3)
        sage: sigma.parent()
        Galois group of Cyclotomic Field of order 8 and degree 4

    Converting from a galois homomorphism to the corresponding
    integer::

        sage: L = CyclotomicField(30)
        sage: sigma = L.galois_group()[2]; sigma
        (1,3)(2,4)(5,7)(6,8)
        sage: cyclotomic_galois_isomorphism(sigma)
        19

    If the modulus is given can convert from a galois homomorphism
    over any field::

        sage: K = QuadraticField(5)
        sage: sigma = K.galois_group().gen()
        sage: cyclotomic_galois_isomorphism(sigma, N=5)
        2
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
