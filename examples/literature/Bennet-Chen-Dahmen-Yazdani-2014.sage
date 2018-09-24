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

E1 = FreyCurve([0, 0, 0, -3*a, -2*c], condition=C1)

K.<sqrt3> = QuadraticField(3)
G = K.galois_group()
E2 = FreyQcurve([0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0],
                isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                condition=C)

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

K.<sqrt3> = QuadraticField(3)
G = K.galois_group()
E21 = FreyQcurve([0, 4*(sqrt3 - 1)*t, 0, -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 + sqrt3)*t^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C21)
E22 = FreyQcurve([0, 4*(sqrt3 - 1)*t, 0, -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 - sqrt3)*t^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C22)
E31 = FreyQcurve([0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 + 3)*s^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C3)
E32 = FreyQcurve([0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 - 3)*s^2), 0],
                 isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                 condition=C3)

