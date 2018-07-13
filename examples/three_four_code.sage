load('/home/jln224/Documents/SageFiles/tmp1.sage')

### Initializing elliptic curves
DA = DiophantineAnalyzer(['a','b'])
DA.add_restriction(CoprimeRestriction([a,b]))
f = (a-b)^4 + a^4 + (a+b)^4
DA.add_restriction(PolynomialPowerRestriction(f,5))
K.<w> = QuadraticField(30)
E1 = EllipticCurve([0, -40*b, 0, -20*(w*a^2 - 10*b^2 + 2*w*b^2), 0])
#E1 = twist_elliptic_curve(E1, w)
E2 = EllipticCurve([0, -60*a, 0, -30*(-15*a^2 + 3*w*a^2 + w*b^2), 0])
E2_orig = E2
E2 = twist_elliptic_curve(E2, 6 + w)
E2t = twist_elliptic_curve(E2, 2*w - 11)
# E2 = twist_elliptic_curve(E2, 9*w - 50)
# E2 = twist_elliptic_curve(E2, w)

### Computing conductors
answer1, answer1_models = DA.compute_conductor(E1, model=True)
answer2 = DA.compute_conductor(E2)

### Model printing:
def print_models(models):
    for p in models:
        print p
        for x in models[p]:
            print x[0]
            print x[1].give_as_congruence_condition()
            print ""
        print ""
    
### Useful function
def action_on_base(f, poly):
    if poly in f.domain():
        return f(poly)
    result = 0
    for mon in poly.monomials():
        result += action_on_base(f,poly.monomial_coefficient(mon)) * mon
    return result


### Isogeny of twist for E1
S = K[a,b]
a2 = S(E1.a2())
a4 = S(E1.a4())
sigma = K.automorphisms()[1]
r = a2^2 - 4 * a4
E1_isog = EllipticCurve([0,-2*a2,0,r,0])
E1_conjug = EllipticCurve([0,action_on_base(sigma,a2),0,action_on_base(sigma,a4),0])

### Special Conductor Computation option 1 curve 1:
K1_d = K
G = CyclotomicField(15).galois_group()
g0, g1 = G.gens()
g = g0^2 * g1
K1_eps = G.subgroup([g^k for k in range(g.order())]).fixed_field()[0]
K1_beta = K1_d.composite_fields(K1_eps)[0]
S1 = K1_beta[a,b]
E1_beta = EllipticCurve([S1(a) for a in E1.a_invariants()])
answer, answer_models = DA.compute_conductor(E1_beta, model=True)

### Special Conductor Computation option 2 curve 1:
E1t = twist_elliptic_curve(E1,w)
K1_d = K
G = CyclotomicField(20).galois_group()
g0, g1 = G.gens()
g = g0^2 * g1
K1_eps = G.subgroup([g^k for k in range(g.order())]).fixed_field()[0]
K1_beta = K1_d.composite_fields(K1_eps)[0]
S1 = K1_beta[a,b]
E1_beta = EllipticCurve([S1(a) for a in E1t.a_invariants()])
answer, answer_models = DA.compute_conductor(E1_beta, model=True)

### Representing the norm of this conductor in a decent way:
delta1 = K1_beta.discriminant()
ramif = [p[0] for p in delta1.factor()]
P = {}
for p in ramif:
    P[p] = K1_beta.primes_above(p)
result = delta1^2
for p in ramif:
    for i in range(len(P[p])):
        prime = P[p][i]
        exp = answer._dict[prime][0][0]
        result = result * prime.absolute_norm()^exp
result = result^(1/4)
print result.factor()

### Isogeny of twist for E2
S = K[a,b]
a2 = S(E2.a2())
a4 = S(E2.a4())
sigma = K.automorphisms()[1]
r = a2^2 - 4 * a4
E2_isog = EllipticCurve([0,-2*a2,0,r,0])
E2_conjug = EllipticCurve([0,action_on_base(sigma,a2),0,action_on_base(sigma,a4),0])

### Special Conductor Computation option 1 curve 2:
K2_d = K
G = CyclotomicField(15).galois_group()
g0, g1 = G.gens()
g = g0^2 * g1
K2_eps = G.subgroup([g^k for k in range(g.order())]).fixed_field()[0]
K2_beta = K2_d.composite_fields(K2_eps)[0]
S2 = K2_beta[a,b]
E2_beta = EllipticCurve([S2(a) for a in E2.a_invariants()])
answer, answer_models = DA.compute_conductor(E2_beta, model=True)

### Special Conductor Computation option 2 curve 2:
E2t = twist_elliptic_curve(E2,w)
K2_d = K
G = CyclotomicField(15).galois_group()
g0, g1 = G.gens()
g = g0^2 * g1
K2_eps = G.subgroup([g^k for k in range(g.order())]).fixed_field()[0]
K2_beta = K2_d.composite_fields(K2_eps)[0]
S2 = K2_beta[a,b]
E2_beta = EllipticCurve([S2(a) for a in E2.a_invariants()])
answer, answer_models = DA.compute_conductor(E2_beta, model=True)

### Representing the norm of this conductor in a decent way:
delta2 = K2_beta.discriminant()
PK = {}
for P in answer._dict.iterkeys():
    p = P.smallest_integer()
    if p not in PK:
        PK[p] = []
    if P not in PK[p]:
        PK[p].append(P)
for p in delta2.prime_factors():
    if p not in PK:
        PK[p] = []
    for P in K2_beta.primes_above(p):
        if P not in PK[p]:
            PK[p].append(P)
result = delta2^2
for p in PK.iterkeys():
    for i in range(len(PK[p])):
        prime = PK[p][i]
        exp = min([tmp[0] for tmp in answer._dict[prime]])
        print p, i, exp
        print len(answer._dict[prime])
        result = result * prime.absolute_norm()^exp
result = result^(1/4)
print result.factor()

### Computing number of points mod p
L = CyclotomicField(8).composite_fields(QuadraticField(2))[0]
R.<T> = L[]
p = 11
N = 4 * p
result = []
G = L.galois_group()
for aa in range(1,N,4):
    for bb in range(1,N,2):
        if not p.divides(gcd(aa,bb)):
            ainv = [an(aa,bb) for an in E2_beta.a_invariants()]
            E_test = EllipticCurve(ainv)
            f = 1
            for prime in K2_beta.primes_above(p):
                Np = prime.norm()
                fp = Np.factor()[0][1]
                ap = 1 + Np - E_test.reduction(prime).count_points()
                f = f * (T^(2*fp) - ap * T^fp + Np)
            dic = dict()
            for factor, exponent in f.factor():
                dic[factor] = exponent
            remaining = []
            while sum(dic.values()) > 0:
                for factor, exponent in dic.iteritems():
                    if exponent > 0:
                        break
                for sigma in G:
                    twisted_factor = R([sigma(coeff) for coeff in factor.list()])
                    dic[twisted_factor] -= 1
                if factor.degree() == 2:
                    ap = -factor.list()[1]
                    flag = False
                    for sigma in G:
                        if sigma(ap) in result:
                            flag = True
                            break
                    if not flag:
                        result.append(ap)
                elif factor.degree() == 1:
                    remaining.append(factor)
                else:
                    print "Error: factor", factor, "is too large."
            for i in range(len(remaining)):
                for j in range(i+1, len(remaining)):
                    factor = remaining[i] * remaining[j]
                    print factor.list()[0]
                    ap = -factor.list()[1]
                    flag = False
                    for sigma in G:
                        if sigma(ap) in result:
                            flag = True
                            break
                    if not flag:
                        result.append(ap)
for ap in result:
    ap_list = []
    for sigma in G:
        if sigma(ap) not in ap_list:
            ap_list.append(sigma(ap))
    print ap_list

### Retrieve labda
def get_lambda(E, K):
    """"""
    a2 = E.a2()
    a4 = E.a4()
    sigma = K.automorphisms()[1]
    r = a2^2 - 4 * a4
    E_isog = EllipticCurve([0,-2*a2,0,r,0])
    E_conjug = EllipticCurve([0,action_on_base(sigma,a2),0,action_on_base(sigma,a4),0])
    gamma = E_isog.a2() / E_conjug.a2()
    if gamma^2 != E_isog.a4() / E_conjug.a4():
        print "Error!"
    gamma = gamma.numerator().constant_coefficient() / gamma.denominator().constant_coefficient()
    return -2 / sqrt(gamma)
