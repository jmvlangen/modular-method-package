def galois_field_extend(sigma, K, embedding=None):
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
    """Changes a Galois action to one on a specified field"""
    L = sigma.parent().number_field()
    M, L_to_M, K_to_M = composite_field(L, K, give_maps=True)
    return galois_field_restrict(galois_field_extend(sigma, M, embedding=L_to_M), K, embedding=K_to_M)

@cached_function
def cyclotomic_galois_isomorphism(s, N=None):
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
