### parameter 2 for a 2-isogeny
E = Qcurve_with_2_isogeny(2)
print "Curve:", E
print "Decomposition field:", E.decomposition_field() # Small field!
Lbeta = E.splitting_image_field()
print "Coefficient field of associated newform:", Lbeta # Field of the same size!
print "Splitting character:", E.splitting_character() # Character is trivial
N = ZZ(sqrt(E.conductor_restriction_of_scalars()))
print "level of newform:", N.factor()

cfs = magma.CuspForms(N)
nfs = magma.Newforms(cfs)
Em = magma(E)
LE = Em.LSeries()

nfs0 = []
for f in nfs:
    Kf = f[1].BaseField().sage()
    if Kf.degree() == 2 and Kf.discriminant().squarefree_part() == 2:
        nfs0.append(f[1])
print len(nfs0)

candidates = copy(nfs0)
print len(candidates)
for p in prime_range(5, 100):
    remove = []
    PE = LE.EulerFactor(p)
    print PE
    for f in candidates:
        Pf = Euler_factor_modular_form(f, p)
        print f.BaseField()
        print Pf, f.Coefficient(p)
        if PE != Pf:
            remove.append(f)
    for f in remove:
        candidates.remove(f)
    if len(candidates) <= 1:
        break
    print ""
print candidates
