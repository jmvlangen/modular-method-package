### Initialization
K.<a> = QuadraticField(5).composite_fields(QuadraticField(-19))[0]
C = K.class_group()
G = K.galois_group()
I = (C.gens()[0]^2)
gamma = ((I.ideal())^2).gens_reduced()[0]
a = {s : sqrt(s(gamma)/gamma) for s in G}
_c = {(s, t) : a[s] * s(a[t]) / a[s*t] for s in G for t in G}
def c(s,t):
    return _c[(s,t)]
U = K.unit_group()
N = min(u.order() for u in U.gens())
Gunit = {s : matrix([U(s(K(U.gens()[i]))).list() for i in range(len(U.gens()))]) for s in G}
