# On the equation a^3 + b^(3*n) = c^2
# Acta Arith. 163 (2014), no. 4, 327--343
R.<s,t> = QQ[]

a = 2*(s^4 + 2*t*s^3 + 2*t^3*s + t^4)
bn = (s-t)^4 - 12*(s*t)^2
c = 3*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)
C1 = (CoprimeCondition(['s','t']) &
      ~CongruenceCondition(s-t, 2) &
      ~CongruenceCondition(s-t, 3) &
      PowerCondition(bn, 5))

E1 = FreyCurve([0, 0, 0, -3*a, -2*c], condition=C1)
K.<sqrt3> = QuadraticField(3)
G = K.galois_group()
E2 = FreyQcurve([0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0],
                isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)},
                condition=C1)
