# On the equation a^3 + b^(3*n) = c^2
# Acta Arith. 163 (2014), no. 4, 327--343
K.<sqrt3> = QuadraticField(3)
R.<s,t> = K[]
bn = (s-t)^4 - 12*(s*t)^2
C = (CoprimeCondition(['s','t']) &
     ~CongruenceCondition(s-t, 2) &
     ~CongruenceCondition(s-t, 3) &
     CongruenceCondition(bn - 1, 2) &
     PowerCondition(bn, 5))
G = K.galois_group()
E2 = FreyQcurve([0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0],
                isogenies={G[0]: (QQ(1), 1), G[1]: (-1 - sqrt3, 2)})
