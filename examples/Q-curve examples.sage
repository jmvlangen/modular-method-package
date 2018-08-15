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

### Wait

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
    
### Tests
for i in range(len(nfs)):
    Kf = nfs[i+1][1].BaseField().sage()
    if Kf.degree() == 4:
        print i
        print Kf.degree()
        print [tmp[0].discriminant().squarefree_part() for tmp in Kf.subfields(degree=2)]
        print ""
