def dirichlet_fixed_field(eps):
    r"""
    Gives the fixed field of the kernel of a Dirichlet chracter.

    INPUT:
    
    - ``eps`` -- A Dirichlet character.

    OUTPUT:
    
    A subfield K of $\Q(\zeta_N)$ such that
     - N is the conductor of eps
     - For each n in the kernel of eps, K is invariant
       under the map $ \zeta_N \mapsto \zeta_N^n $.
    """
    eps = eps.primitive_character()
    N = eps.conductor()
    return fixed_field([cyclotomic_galois_isomorphism(s, N) for s in eps.kernel()])
    # H = [cyclotomic_galois_isomorphism(s, N) for s in eps.kernel()]
    # G = CyclotomicField(N).galois_group()
    # H = G.subgroup(H)
    # if H == G:
    #     return QQ
    # result = H.fixed_field()
    # if isinstance(result, tuple):
    #     return result[0]
    # else:
    #     return result

def dirichlet_to_galois(eps):
    r"""
    Gives the galois character associated to a Dirichlet character.

    INPUT:
    
    - ``eps`` -- A Dirichlet character.

    OUTPUT:
    
    A function f that is the galois version of the Dirichlet character
    eps, e.g. we have
    ..MATH:
    
        f( \zeta_N \mapsto \zeta_N^n ) = eps(n)

    where $N$ is any number divisible by the conductor of eps. Note that
    this function also gives output on any element that is not necessary
    of the right galois group, by taking an extension and restriction
    via the common composite field.
    """
    eps = eps.primitive_character()
    N = eps.conductor()
    def eps_galois(sigma):
        return eps(cyclotomic_galois_isomorphism(sigma, N=N))
    return eps_galois

def character_for_root(a):
    r"""
    Gives the Kronecker character of a given number.

    INPUT:

    - A non-zero integer a

    OUTPUT:

    A Dirichlet character eps such that the associated galois
    character has fixed field $ \Q(\sqrt{a}) $
    """
    a = a.squarefree_part()
    ls = a.prime_factors()
    ls.append(a.sign())
    N = 1
    eps = DirichletGroup(1)[0]
    for i in range(len(ls)-1):
        m = mod(ls[i],4)
        if m != 2:
            D = DirichletGroup(ls[i])
        if m == 1:
            N *= ls[i]
            eps = eps.extend(N) * (D.gen()^ZZ(D.order()/2)).extend(N)
        elif m == 2:
            N *= 8
            eps = eps.extend(N) * _character_for_2(2).extend(N)
        elif m == 3:
            N *= ls[i]
            ls[-1] *= -1
            eps = eps.extend(N) * (D.gen()^ZZ(D.order()/2)).extend(N)
        else:
            raise Exception("Invalid prime factor %s"%ls[i])
    if ls[-1] == -1:
        if not 4.divides(N):
            N *= 4
        D = DirichletGroup(4)
        eps = eps.extend(N) * D.gen().extend(N)
    return eps

def _character_for_2(a):
    r"""
    Gives the Kronecker character of some numbers that have to do with 2.

    INPUT:

    - An integer from the list [-1, 2, -2]

    OUTPUT:

    A Dirichlet character eps such that the associated galois
    character has fixed field $ \Q(\sqrt{a}) $
    """
    D = DirichletGroup(8)
    if a == -1:
        b = 5
    elif a == 2:
        b = 7
    elif a == -2:
        b = 3
    else:
        raise Exception("Value %s not accepted"%a)
    for eps in D:
        if len(eps.kernel()) == 2 and b in eps.kernel():
            return eps
    raise Exception("No character found!")
