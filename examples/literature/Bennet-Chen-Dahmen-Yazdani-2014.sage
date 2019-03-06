# On the equation a^3 + b^(3*n) = c^2
# Acta Arith. 163 (2014), no. 4, 327--343

### Case c odd
R.<s,t> = QQ[]
a = 2*(s^4 + 2*t*s^3 + 2*t^3*s + t^4)
bn = (s-t)^4 - 12*(s*t)^2
c = 3*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)

C = (CoprimeCondition(['s','t']) &
      ~CongruenceCondition(s-t, 2) &
      ~CongruenceCondition(s-t, 3) &
      PowerCondition(bn, 5))

E1 = FreyCurve([0, 0, 0, -3*a, -2*c], condition=C)
print ""; print E1.conductor()
# 64*3^n0*Rad_P( (-1728) * (s^4 - 4*s^3*t - 6*s^2*t^2 - 4*s*t^3 + t^4)^3 )
# where 
# n0 = 2 if ('s', 't') == (1, 2), (2, 1) mod 3
#      3 if ('s', 't') is 1 of 4 possibilities mod 3

K.<sqrt3> = QuadraticField(3)
G = K.galois_group()
E2 = FreyQcurve([0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0],
                isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                condition=C)
print ""; print E2.conductor()
# (64)*Rad_P( ((-960*sqrt3 + 1664)) * (s^2 + (2*sqrt3 - 2)*s*t + t^2) * (s^2 + (-2*sqrt3 - 2)*s*t + t^2)^2 )
print ""; print E2.conductor_restriction_of_scalars()
# 589824*Norm(Rad_P( ((-960*zeta0 + 1664)) * (s^2 + (2*zeta0 - 2)*s*t + t^2) * (s^2 + (-2*zeta0 - 2)*s*t + t^2)^2 ))
print ""; print E2.newform_levels()
# [(768,)]

### Case c even:
R.<s,t> = QQ[]
a1 = 3*s^4 + 6*t^2*s^2 - t^4
a2 = -3*s^4 + 6*t^2*s^2 + t^4
bn11 = -3*s^4 + 6*t^2*s^2 + t^4
bn12 = 3*s^4 + 6*t^2*s^2 - t^4
bn21 = 3*(s^2 + t^2)^2 - 4*t^4
bn22 = 3*(s^2 - t^2)^2 - 4*t^4
bn3 = (t^2 + 3*s^2)^2 - 12*s^4
c = 6*s*t*(3*s^4 + t^4)

C = (CoprimeCondition(['s', 't']) &
     ~CongruenceCondition(t, 3) &
     ~CongruenceCondition(s-t, 2))
C11 = C & PowerCondition(bn11, 5)
C12 = C & PowerCondition(bn12, 5)
C21 = C & PowerCondition(bn21, 5)
C22 = C & PowerCondition(bn22, 5)
C3 = C & PowerCondition(bn3, 5)

E11 = FreyCurve([0, 0, 0, -3*a1, -2*c], condition=C11)
E12 = FreyCurve([0, 0, 0, -12*a2, -16*c], condition=C12)
print ""; print E11.conductor()
# 32*3^n0*Rad_P( (1728) * (3*s^4 - 6*s^2*t^2 - t^4)^3 )
#  where 
# n0 = 2 if ('s', 't') == (0, 1), (0, 2) mod 3
#      3 if ('s', 't') is 1 of 4 possibilities mod 3
print ""; print E12.conductor()
# 32*3^n0*Rad_P( (-110592) * (3*s^4 + 6*s^2*t^2 - t^4)^3 )
#  where 
# n0 = 2 if ('s', 't') == (0, 1), (0, 2) mod 3
#      3 if ('s', 't') is 1 of 4 possibilities mod 3

K.<sqrt3> = QuadraticField(3)
G = K.galois_group()
E21 = FreyQcurve([0, 4*(sqrt3 - 1)*t, 0, -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 + sqrt3)*t^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C21)
E22 = FreyQcurve([0, 4*(sqrt3 - 1)*t, 0, -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 - sqrt3)*t^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C22)
print ""; print E21.conductor()
# (64)*Rad_P( ((39936*sqrt3 - 69120)) * (s^2 + (2/3*sqrt3 + 1)*t^2) * (s^2 + (-2/3*sqrt3 + 1)*t^2)^2 )
print ""; print E22.conductor()
# (64)*Rad_P( ((39936*sqrt3 - 69120)) * (s^2 + (2/3*sqrt3 - 1)*t^2) * (s^2 + (-2/3*sqrt3 - 1)*t^2)^2 )
print ""; print E21.conductor_restriction_of_scalars()
# 589824*Norm(Rad_P( ((39936*zeta0 - 69120)) * (s^2 + (2/3*zeta0 + 1)*t^2) * (s^2 + (-2/3*zeta0 + 1)*t^2)^2 ))
print ""; print E22.conductor_restriction_of_scalars()
# 589824*Norm(Rad_P( ((39936*zeta0 - 69120)) * (s^2 + (2/3*zeta0 - 1)*t^2) * (s^2 + (-2/3*zeta0 - 1)*t^2)^2 ))
print ""; print E21.newform_levels()
# [(768,)]
print ""; print E22.newform_levels()
# [(768,)]

E31 = FreyQcurve([0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 + 3)*s^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C3)
E32 = FreyQcurve([0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 - 3)*s^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C3)
print ""; print E31.conductor()
# (192)*Rad_P( ((-1492992*sqrt3 + 2612736)) * (s^2 + (-2/3*sqrt3 - 1)*t^2) * (s^2 + (2/3*sqrt3 - 1)*t^2)^2 )
print ""; print E32.conductor()
# (192)*Rad_P( ((-20901888*sqrt3 + 36205056)) * (s^2 + (-2/3*sqrt3 + 1)*t^2) * (s^2 + (2/3*sqrt3 + 1)*t^2)^2 )
print ""; print E31.conductor_restriction_of_scalars()
# 5308416*Norm(Rad_P( ((-1492992*zeta0 + 2612736)) * (s^2 + (-2/3*zeta0 - 1)*t^2) * (s^2 + (2/3*zeta0 - 1)*t^2)^2 ))
print ""; print E32.conductor_restriction_of_scalars()
# 5308416*Norm(Rad_P( ((-20901888*zeta0 + 36205056)) * (s^2 + (-2/3*zeta0 + 1)*t^2) * (s^2 + (2/3*zeta0 + 1)*t^2)^2 ))
print ""; print E31.newform_levels()
# [(2304,)]
print ""; print E32.newform_levels()
# [(2304,)]
