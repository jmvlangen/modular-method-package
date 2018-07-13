load("~/Documents/SageFiles/tmp1.sage")
###Example 1
DA = DiophantineAnalyzer(['u','v'])
f = u^2 + v^2
s = v^2
t = u^2
DA.add_restriction(CoprimeRestriction([u,v]))
DA.add_restriction(PolynomialPowerRestriction(f, 5))
K.<r5> = QuadraticField(-2)
delta = (-5 + 3*r5)/2
a4s = -3 * delta * ((3 + 2*r5)*s - 3*t)
a6s = 4*v*((17 - 4*r5)*s - (45 - 18*r5)*t)
Es = EllipticCurve([0,0,0,a4s,a6s])
a4t = -3*2^2*r5*(3*s - (15 - 10*r5)*t)
a6t = 2^5*5*u*(9*s - (45 - 14*r5)*t)
Et = EllipticCurve([0,0,0,a4t,a6t])
# answer = DA.compute_conductor(E, verbose=True)
