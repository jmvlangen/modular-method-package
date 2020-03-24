def _poly_galois_conjugates(group, poly):
    for sigma in group:
        yield poly.parent()([sigma(coeff) for coeff in poly.list()])

def _galois_orbit_factors(group, poly):
    result = []
    dic = dict()
    for factor, exponent in poly.factor():
        dic[factor] = exponent
    result.append(factor)
    while sum(dic.values()) > 0:
        for factor in dic:
            if dic[factor] > 0:
                break
        for twisted_factor in _poly_galois_conjugates(group, factor):
            dic[twisted_factor] -= 1
    return result

def _char_poly_p_factor(E, prime, Np, fp, T):
    if E.has_good_reduction(prime):
        ap = 1 + Np - E.reduction(prime).count_points()
        c0 = Np
    elif E.has_split_multiplicative_reduction(prime):
        ap = 1
        c0 = 0
    elif E.has_multiplicative_redution(prime):
        ap = -1
        c0 = 0
    else:
        ap = 0
        c0 = 0
    return T^(2*fp) - ap * T^fp + c0

def _char_poly_p(E, primes, Nps, fps, T):
    result = 1
    for i in range(len(primes)):
        result *= _char_poly_p_factor(E, primes[i], Nps[i], fps[i], T)
    return result

def _prime_data(field, p):
    primes = field.primes_above(p)
    Nps = [prime.norm() for prime in primes]
    fps = [Np.factor()[0][1] for Np in Nps]
    return primes, Nps, fps

def _ngens(n, m):
    r = [0]*m
    i = 0
    while i < m:
        if i == 0:
            yield tuple(r)
        r[i] += 1
        if r[i] >= n:
            r[i] = 0
            i += 1
        else:
            i = 0

def _mod_tuple(l, mod):
    result = []
    for n in l:
        result.append(n.mod(mod))
    return tuple(result)

def _elliptic_possibilities(E, p, vals, mod):
    N = p * mod
    n = len(E.base().variable_names())
    for v in _ngens(N, n):
        if _mod_tuple(v, mod) in vals and not p.divides(gcd(v)):
            ainv = [a(v) for a in E.a_invariants()]
            yield EllipticCurve(ainv)

def _add_galois_orbit(ls, number, group):
    flag = True
    for sigma in group:
        if sigma(number) in ls:
            flag = False
            break
    if flag:
        ls.append(number)

def _orbit_list_to_element_list(orbit_list, G):
    element_list = []
    for x in orbit_list:
        for sigma in G:
            if sigma(x) not in element_list:
                element_list.append(sigma(x))
    return element_list

def ap_possibilities(E, p, K, L, vals, mod, c0):
    G = L.galois_group()
    R = PolynomialRing(L,name='T')
    T = R.gens()[0]
    result = []
    primes, Nps, fps = _prime_data(K, p)
    for Ei in _elliptic_possibilities(E, p, vals, mod):
        factors = _galois_orbit_factors(G, _char_poly_p(Ei, primes, Nps, fps, T))
        remaining = []
        for factor in factors:
            if factor.degree() == 2:
                for twisted_factor in _poly_galois_conjugates(G, factor):
                    if twisted_factor.list()[0] == c0 and twisted_factor.list()[1] not in result:
                        result.append(twisted_factor.list()[1])
            elif factor.degree() == 1:
                for twisted_factor in _poly_galois_conjugates(G, factor):
                    remaining.append(twisted_factor)
            else:
                print("Error:", factor, "has an incorrect degree!")
        for i in range(len(remaining)):
            for j in range(i+1, len(remaining)):
                factor = remaining[i] * remaining[j]
                if factor.list()[0] == c0:
                    if factor.list()[1] not in remaining:
                        result.append(factor.list()[1])
    return result

def ap_possibilities_print(E, p_start, p_stop, K, L, vals, mod, eps):
    for p in prime_range(p_start, p_stop):
        print(p, ap_possibilities(E, p, K, L, vals, mod, eps(p)*p))

def eliminate_nf(E, K_E, g, K_f, eps, vals, mod, bad_p, max_p=200, l_cap=[]):
    K_g = g.base_ring()
    L, K_g_to_L, K_f_to_L, k = K_g.composite_fields(K_f, names='b', both_maps=True)[0]
    possible_l = 0
    important_p = []
    for p in prime_range(max_p):
        if p not in bad_p:
            if possible_l == 0:
                prev_pos = -1
            else:
                prev_pos = len(possible_l)
            ap_f_list = ap_possibilities(E, p, K_E, K_f, vals, mod, eps(p)*p)
            print(ap_f_list)
            ap_g = g.coefficient(7)
            nrm = 1
            for ap_f in ap_f_list:
                nrm *= (K_f_to_L(ap_f) - K_g_to_L(ap_g)).absolute_norm()
            if nrm != 0:
                l_list = [f for f,e in nrm.factor()]
                if possible_l == 0:
                    possible_l = l_list
                else:
                    for l in possible_l:
                        if l not in l_list:
                            possible_l.remove(l)
            # Printing a helpful message for the user
            if possible_l != 0 and prev_pos != len(possible_l):
                print(p, possible_l)
                important_p.append(p)
            else:
                print(p, "No change")
            # Stopping if cap is reached, i.e. can not do any better.
            if possible_l != 0:
                flag = True
                for l in possible_l:
                    if l not in l_cap:
                        flag = False
                        break
                if flag:
                    break
    return possible_l, important_p
