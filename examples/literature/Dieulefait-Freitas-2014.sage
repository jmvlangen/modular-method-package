# The fermat-type equations x^5 + y^5 = 2 z^p or 3 z^p solved through Q-curves.
# Mathematics of Computation, volume 83 (2014), number 286, pages 917-933

R.<a,b> = QQ[]
dcp = a^5 + b^5
phi = a + b
i = 0
while phi == a + b:
    phi = dcp.factor()[i][0]
    i += 1

K.<sqrt5> = QuadraticField(5)
(phi2, e2), (phi1, e1o) = phi.change_ring(K).factor()
w = (-1 + sqrt(5))/2
G.<sigma> = K.galois_group()

C = (CoprimeCondition([a,b]))

L.<sqrtm2> = QuadraticField(-2)
E = FreyQcurve([0, 2*(a + b), 0, - sigma(w)*phi1, 0],
               isogenies={G(1): (QQ(1), 1), sigma: (sqrtm2, 2)},
               condition=C)

S.<x> = QQ[]
K.<theta> = NumberField(x^4 - 5*x^2 + 5)
gamma = 2*theta^2 - theta - 5
Eg = E.twist(gamma)

# Field is bigger than necessary
