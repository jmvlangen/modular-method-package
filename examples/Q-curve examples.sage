### parameter 2 for a 2-isogeny
E = Qcurve_with_2_isogeny(2)
print "Original curve:", E
print "Decomposition field:", E.decomposition_field()

print ""

E = E.complete_definition_twist([2])
print "Better choice curve:", E
print "Degree field:", E.degree_field()
print "Decomposition_field:", E.decomposition_field() # Is Kd
answer = E.does_decompose()
print "Does it decompose?", answer

print ""

E = E.decomposable_twist()
print "Best choice:", E
print "Degree field:", E.degree_field()
print "Decomposition_field:", E.decomposition_field() # Is Kd
answer = E.does_decompose()
print "Does it decompose?", answer

Lbeta = E.splitting_image_field()
print "newform coefficient field:", Lbeta # Of degree 2!
N = ZZ(sqrt(E.conductor_restriction_of_scalars()))
print "level of newform:", N.factor()
print "character:", E.splitting_character()^(-1)

cfs = magma.CuspForms(magma(N))
nfs = magma.Newforms(cfs)
Em = magma(E)

nfs0 = []
for f in nfs:
    Kf = f[1].BaseField().sage()
    if Kf.degree() == Lbeta.degree() and \
       Kf.discriminant().squarefree_part() == Lbeta.discriminant().squarefree_part():
        nfs0.append(f[1])

candidates = copy(nfs0)
for p in prime_range(5, 100):
    remove = []
    PE = Em.EulerFactor(p)
    for i in range(len(candidates)):
        f = candidates[i]
        Pf = Euler_factor_modular_form(f, p)
        if PE != Pf:
            remove.append(f)
    for f in remove:
        candidates.remove(f)
    if len(candidates) <= 1:
        break

if len(candidates) >= 1:
    f = candidates[0]
    print "Associated newform:", f
    print "Verified for primes up to 1000:", all(Em.EulerFactor(p) == Euler_factor_modular_form(f, p) for p in prime_range(5, 1000))
else:
    print "Error: No corresponding newform found"

# K = E.decomposition_field()
# f = candidates[0]
# for p in prime_range(5, 100):
#     print p, len(K.primes_above(p)), p % 8
#     print Em.EulerFactor(p)
#     print Euler_factor_modular_form(f, p)
#     print Em.EulerFactor(p) == Euler_factor_modular_form(f, p)
#     print ""

### Tests
R.<x> = QQ[]
for p in prime_range(5,100):
    a0, a1, a2, a3, a4 = Em.EulerFactor(p).Coefficients().sage()
    print p
    if a1 != 0:
        epsp = a3 / (a1*p)
        print "eps(p) =", epsp
        print "ap minpoly =", R([a2 - 2*epsp*p, a1, 1]).factor()
    else:
        print "Two cases:"
        epsp = 1
        print "eps(p) =", epsp
        print "ap minpoly =", R([a2 - 2*epsp*p, a1, 1]).factor()
        epsp = -1
        print "eps(p) =", epsp
        print "ap minpoly =", R([a2 - 2*epsp*p, a1, 1]).factor()
    print ""

### 2
ls = []
for p in prime_range(5,1000):
    a0, a1, a2, a3, a4 = Em.EulerFactor(p).Coefficients().sage()
    if a1 == 0 and (a2 + 2*p).squarefree_part() == 2:
        ls.append(p)

for n in range(1, 20):
    mod_ls = []
    for p in ls:
        if p % n not in mod_ls:
            mod_ls.append(p % n)
    flag = True
    for p in prime_range(5,1000):
        if p % n in mod_ls and p not in ls:
            flag = False
            break
    if flag:
        print n
