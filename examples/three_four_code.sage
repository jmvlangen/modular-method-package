K.<w> = QuadraticField(30)
R.<a,b> = ZZ[]
f = (a-b)^4 + a^4 + (a+b)^4
C = CoprimeCondition([a,b]) & PowerCondition(f, 5)
G = K.galois_group()

### The curves
Lm2.<sqrtm2> = QuadraticField(-2)
E1 = FreyQcurve([0, -40*b, 0, -20*(w*a^2 - 10*b^2 + 2*w*b^2),0],
                condition=C,
                isogenies={G[0]: (QQ(1),1), G[1]: (sqrtm2,2)})
E1 = E1.decomposable_twist()
E2 = FreyQcurve([0, -60*a, 0, -30*(-15*a^2 + 3*w*a^2 + w*b^2),0],
                condition=C,
                isogenies={G[0]: (QQ(1),1), G[1]: (sqrtm2,2)})
E2 = E2.decomposable_twist()

### the data
print ""; print E1._newform_levels()
# [(15360, 15360, 76800, 76800), (76800, 76800, 15360, 15360)]
print ""; print E2._newform_levels()
# [(23040, 23040, 115200, 115200), (115200, 115200, 23040, 23040)] if ('a', 'b') == (1, 0) mod 2 and ('a', 'b') == (1, 0) mod 2
# []                                                               if ('a', 'b') == (1, 0) mod 2 and ('a', 'b') == (1, 1) mod 2
# []                                                               if ('a', 'b') == (1, 1) mod 2 and ('a', 'b') == (1, 0) mod 2
# [(11520, 11520, 57600, 57600), (57600, 57600, 11520, 11520)]     if ('a', 'b') == (1, 1) mod 2 and ('a', 'b') == (1, 1) mod 2

### Computing newspaces
eps = E1.splitting_character()
# eps = E2.splitting_character() # Note that they are in fact the same!!!
Dm = magma.DirichletGroup(eps.conductor(), magma(eps.base_ring()))
for eps_m in Dm.Elements():
    if (eps_m(11) == eps(11) and eps_m(7) == eps(7)): # 11 -> -1 , 7 -> zeta4, conductor 15
        break

cfs1 = magma.CuspForms(magma.DirichletGroup(15360, magma(eps.base_ring()))(eps_m))
cfs21 = magma.CuspForms(magma.DirichletGroup(11520, magma(eps.base_ring()))(eps_m))
cfs22 = magma.CuspForms(magma.DirichletGroup(23040, magma(eps.base_ring()))(eps_m))

