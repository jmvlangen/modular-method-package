# Solving Fermat-type equations via modular Q-curves over polyquadratic fields
# ???

R.<a,b> = QQ[]
cp = a^4 + 2*b^2

C = (CoprimeCondition([a,b]) &
     PowerCondition(cp, 5))

K.<r> = QuadraticField(-2)
G.<sigma> = K.galois_group()
E = FreyQcurve([0, 4*a, 0, 2*(a^2 + r*b), 0],
               isogenies={G(1): (QQ(1), 1), sigma: (r, 2)},
               condition=C)
