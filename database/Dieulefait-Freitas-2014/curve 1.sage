# Article: The fermat-type equations x^5 + y^5 = 2 z^p or 3 z^p solved through Q-curves.
# Authors: Luis Dieulefait and Nuno Freitas
# Journal: Mathematics of Computation, volume 83 (2014), number 286, pages 917-933
#
# Equation: a^5 + b^5 = 2*c^l
# Assumptions:
#   - a, b and 2*c are coprime

# Variables
R.<a,b> = QQ[]
dcp = a^5 + b^5
phi = R(dcp / (a + b))
S.<x> = QQ[]
K.<w> = NumberField(x^2 + x - 1)
G.<sigma> = K.galois_group()
phi1 = a^2 + w*a*b + b^2
phi2 = a^2 + sigma(w)*a*b + b^2

# Conditions
C = CoprimeCondition([a,b])

# The curve
# Y^2 = X^3 + 2*(a + b)*X^2 - sigma(w)*phi1*X
a_invariants = [0, 2*(a + b), 0, -sigma(w)*phi1, 0]
L.<sqrtm2> = QuadraticField(-2)
isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm2, 2)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               condition=C)

# Twisting
E = E.decomposable_twist()

# Conductor
print ""
print E.conductor()
# (2, 1/16*wsqrtm2zeta00^3 + 1/4*wsqrtm2zeta00^2 - 3/2*wsqrtm2zeta00 - 4)^n0*(5, 1/8*wsqrtm2zeta00^2 - 2)^n1*Rad_P( ((-5540*wsqrtm2zeta00^3 - 3960*wsqrtm2zeta00^2 + 136800*wsqrtm2zeta00 + 120320)) * (a^2 + (-1/16*wsqrtm2zeta00^3 - 1/8*wsqrtm2zeta00^2 + 3/2*wsqrtm2zeta00 + 2)*a*b + b^2) * (a^2 + (1/16*wsqrtm2zeta00^3 + 1/8*wsqrtm2zeta00^2 - 3/2*wsqrtm2zeta00 - 3)*a*b + b^2)^2 )
#  where 
# n0 =  8 if ('a', 'b') is 1 of 6 possibilities mod 4
#       6 if ('a', 'b') is 1 of 4 possibilities mod 4
#       0 if ('a', 'b') is 1 of 4 possibilities mod 8
#       4 if ('a', 'b') is 1 of 4 possibilities mod 8
# n1 =  2 if ('a', 'b') is 1 of 20 possibilities mod 5
#       0 if ('a', 'b') is 1 of 4 possibilities mod 5
# (article agrees at the prime above 5)
# (article has completely different case for the prime above 2)

# Newforms
print ""
print E.newforms(algorithm='magma')
#
# (article: None should remain, but different elimination methods)
