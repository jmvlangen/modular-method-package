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

### parameter 6 for a 2-isogeny
E = Qcurve_with_2_isogeny(6)
print "Original curve:", E
print "Decomposition field:", E.decomposition_field()

print ""

E = E.complete_definition_twist([2])
print "Better choice curve:", E
print "Degree field:", E.degree_field()
print "Decomposition_field:", E.decomposition_field() # Is not Kd!
answer = E.does_decompose()
print "Does it decompose?", answer

print ""

E = E.decomposable_twist()
print "Best choice:", E
print "Degree field:", E.degree_field()
print "Decomposition_field:", E.decomposition_field() # Is not Kd
answer = E.does_decompose()
print "Does it decompose?", answer

Lbeta = E.splitting_image_field()
print "newform coefficient field:", Lbeta # Of degree 4 as is the decomposition field!
N = ZZ(E.conductor_restriction_of_scalars().nth_root(4))
print "level of newform:", N.factor()
eps = E.splitting_character()^(-1)
Dm = magma.DirichletGroup(eps.conductor())
for eps_m in Dm.Elements():
    flag = True
    for i in range(1, eps.conductor()+1):
        i = ZZ(i)
        flag = (eps(i) == eps_m(i))
        if not flag: # Not the right one
            break
    if flag: # Found the right one
        break
eps_m = magma.DirichletGroup(N)(eps_m)
print "character:", eps_m

try:
    f = nfs[1][1]
    if f.parent() != magma or f.DirichletCharacter() != eps_m:
        raise TypeError("Not the right object.")
except (NameError, IndexError, TypeError):
    cfs = magma.CuspForms(eps_m)
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

### parameter -6 for a 2-isogeny
E = Qcurve_with_2_isogeny(-6)
print "Original curve:", E
print "Decomposition field:", E.decomposition_field()
answer = E.does_decompose()
print "Does it decompose?", answer

print ""

E = E.complete_definition_twist([2])
print "Better choice curve:", E
print "Degree field:", E.degree_field()
print "Decomposition_field:", E.decomposition_field() # Is not Kd!
answer = E.does_decompose()
print "Does it decompose?", answer

print ""

E = E.decomposable_twist()
print "Best choice:", E
print "Degree field:", E.degree_field()
print "Decomposition_field:", E.decomposition_field() # Is not Kd
answer = E.does_decompose()
print "Does it decompose?", answer

Lbeta_ls = E.splitting_image_field('conjugacy')
print "newform coefficient fields:", Lbeta_ls # Two times, of degree 4!
print "Possible levels of the newform:"
N_ls = E._newform_levels()
for N in N_ls:
    print N

### Let's make some attempts at doing this!
twists0 = E.twist_character('conjugacy')
M = lcm(chi.modulus() for chi in twists0)
twists0 = [chi.extend(M) for chi in twists0]
for N in N_ls:
    i_min = min(enumerate(N), key=(lambda x : x[1]))[0] # newform with smallest level
    print "level of newform:", N[i_min].factor()
    chi = twists0[i_min]
    twists = [chi_j * chi^(-1) for chi_j in twists0]
    eps = E.splitting_character('conjugacy')[i_min]
    Dm = magma.DirichletGroup(eps.conductor())
    for eps_m in Dm.Elements():
        candidate = True
        for i in range(1, eps.conductor()+1):
            i = ZZ(i)
            candidate = (eps(i) == eps_m(i))
            if not candidate: # Not the right one
                break
        if candidate: # Found the right one
            break
    eps_m = magma.DirichletGroup(N[i_min])(eps_m)
    print "character:", eps_m
    Lbeta = Lbeta_ls[i_min]
    print "Coefficient field:", Lbeta

    Em = magma(E)
    try:
        f = nfs[1][1]
        if f.parent() != magma or f.DirichletCharacter() != eps_m:
            raise TypeError("Not the right object.")
    except (NameError, IndexError, TypeError):
        cfs = magma.CuspForms(eps_m)
        nfs = magma.Newforms(cfs)
        print "Newforms calculated!"

    nfs0 = []
    for f in nfs:
        Kf = f[1].BaseField().sage()
        if Kf.degree() == Lbeta.degree() and \
           Kf.discriminant().squarefree_part() == Lbeta.discriminant().squarefree_part():
            nfs0.append(f[1])

    candidates = copy(nfs0)
    for p in prime_range(1, 100):
        if not p.divides(N[i_min]):
            remove = []
            PE = Em.EulerFactor(p)
            for f in candidates:
                Pf = Euler_factor_modular_form(f, p, twists=twists)
                if PE != Pf:
                    remove.append(f)
            for f in remove:
                candidates.remove(f)
            if len(candidates) <= 1:
                break

    if len(candidates) >= 1 and all(Em.EulerFactor(p) == Euler_factor_modular_form(candidates[0], p, twists=twists) for p in prime_range(1000) if not p.divides(N[i_min])):
        f = candidates[0]
        print "Associated newform:", f
        print "Verified for primes up to 1000!"
        break
    else:
        print "No corresponding newform found at this level"
    
### Tests
print roots_in_field(Lbeta)
for i in range(len(nfs)):
    Kf = nfs[i+1][1].BaseField().sage()
    if Kf.degree() == 4:
        print i
        print Kf.degree()
        print roots_in_field(Kf)
        print ""

### 2
for p in prime_range(1, 100):
    if True: #not p.divides(N[i_min]):
        print p
        PE = Em.EulerFactor(p)
        print PE
        for f in nfs0:
            Pf = Euler_factor_modular_form(f, p, twists=twists)
            if PE == Pf:
                print Pf, "*"
            else:
                print Pf
        print ""

