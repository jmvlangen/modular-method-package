load("~/Documents/SageFiles/tmp1.sage")
###Example 1
DA = DiophantineAnalyzer(['A','B'])
f = A^4 + 2 * B^2
DA.add_restriction(CoprimeRestriction([A,B]))
DA.add_restriction(PolynomialPowerRestriction(f, 5))
K.<r> = QuadraticField(-2)
E = EllipticCurve([0,4*A,0,2*(A^2 + r * B),0])
answer = DA.compute_conductor(E)
# r^(e1) * Rad_P((A^2 + (-r)*B) * (A^2 + (r)*B)^2)
# where P = {r}
# with e1 = 12 if (A, B) = (1, 1) (mod 2)
#           10 if (A, B) = (1, 0) (mod 2)

###Example 2
DA = DiophantineAnalyzer(['A','B'])
f = A^4 + 3 * B^2
DA.add_restriction(CoprimeRestriction([A,B]))
DA.add_restriction(PolynomialPowerRestriction(f, 5))
K.<r> = QuadraticField(-2)
E = EllipticCurve([0,4*A,0,2*(A^2 + r * B),0])
answer = DA.compute_conductor(E)
# r^(e1) * Rad_P((A^2 + (-r)*B) * (A^2 + (r)*B)^2)
# where P = {r}
# with e1 = 12 if (A, B) = (0, 1) (mod 2)
#           10 if (A, B) = (1, 0) (mod 2)
R.<x> = K[]
f = x^2 - 6
L.<gamma> = K.extension(f)
r = L(r)
gamma = L(gamma)
L = L.absolute_field(names='a')
r = L(r)
gamma = L(gamma)
# Unclear how to twist
