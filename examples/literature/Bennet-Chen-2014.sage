# Multi-Frey Q-curves and the Diophantine equation a^2 + b^6 = c^n
# Algebra & Number Theory volume 6 (2012), no.4

R.<a,b> = QQ[]
cp = a^2 + b^6

C = (CoprimeCondition([a,b]) &
     PowerCondition(cp, 5))

K.<i> = QuadraticField(-1)
G = K.galois_group()
E = FreyQcurve([0, 0, 0, -3*(5*b^3 + 4*a*i)*b, 2*(11*b^6 + 14*i*b^3*a - 2*a^2)],
               isogenies={G[0]: (QQ(1), 1), G[1]: (i , 2)}, #TODO!!!!
               condition=C)