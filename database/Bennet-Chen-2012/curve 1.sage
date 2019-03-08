# Article: Multi-Frey Q-curves and the Diophantine equation a^2 + b^6 = c^n
# Authors: Michael A. Bennett and Imin Chen
# Journal: Algebra & Number Theorey, volume 6 (2012), no. 4
# Source: http://people.math.sfu.ca/~ichen/pub/BeCh2.pdf
#
# Equation: a^2 + b^2 = c^l
# Assumptions:
#   - a, b and c are coprime
#   - l >= 3

# Variables
R.<a,b> = QQ[]
cl = a^2 + b^6 # c^l

# Conditions
C = (CoprimeCondition([a,b]) &
     PowerCondition(cl, 3))

# The curve
K.<i> = QuadraticField(-1)
L.<sqrtm3> = QuadraticField(-3)
G.<sigma> = K.galois_group()
a_invariants = [0, 0, 0, -3*(5*b^3 + 4*a*i)*b, 2*(11*b^6 + 14*i*b^3*a - 2*a^2)]
isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm3 , 3)}
E = FreyQcurve(a_invariants, isogenies=isogenies, condition=C)

# Checking properties described in the article
print (E.j_invariant() ==
       432*i*(b^3*(4*a - 5*i*b^3)^3)/((a - i*b^3)*(a + i*b^3)^3))
print E.discriminant() == -2^8*3^3*(a - i*b^3)*(a + i*b^3)^3
print E.degree_map_image() == [1, 3]
print E.degree_field().is_isomorphic(QQ[sqrt(-1)])
print E.dual_basis() == ([-1], [3])
print (E.splitting_character() ==
       DirichletGroup(4).gen().extend(12) * DirichletGroup(3).gen().extend(12))
print E.splitting_character_field().is_isomorphic(QQ[sqrt(3)])
print E.splitting_field().is_isomorphic(QQ[sqrt(3), sqrt(-1)])
print E.splitting_image_field().is_isomorphic(QQ[sqrt(3), sqrt(-1)])

# Twisting
Kb = E.decomposition_field()
gamma = (-3 + sqrt(Kb(-3))) / 2
Eb = E.twist(gamma)
print Eb.does_decompose()
print Eb.j_invariant() == E.decomposable_twist().j_invariant()

# Conductor
q2 = Kb.prime_above(2)
q3 = Kb.prime_above(3)
N = Eb.conductor(additive_primes=[q2, q3])
print N
# (4)*(1/4*izeta0^2 + 1)^n0*Rad_P( (-186624) * (b^3 + (-1/8*izeta0^3)*a) * (b^3 + (1/8*izeta0^3)*a)^3 )
#  where 
# n0 = 0 if ('a', 'b') is 1 of 24 possibilities mod 9
#      4 if ('a', 'b') is 1 of 48 possibilities mod 9
# (agrees with article)
NR = Eb.conductor_restriction_of_scalars(additive_primes=[q2, q3])
print NR
# 65536*3^(2*n0+4)*Norm(Rad_P( (-186624) * (b^3 + (-1/28*izeta00zeta0^3 - 9/28*izeta00zeta0)*a) * (b^3 + (1/28*izeta00zeta0^3 + 9/28*izeta00zeta0)*a)^3 ))
#  where 
# n0 = 0 if ('a', 'b') is 1 of 24 possibilities mod 9
#      4 if ('a', 'b') is 1 of 48 possibilities mod 9

# Newforms
Eb.number_of_conjugacy_classes()
print ""
nfs = E.newforms(algorithm='magma')
print nfs
#
# (article: the only form at level 48 is CM and 2 forms at level 432 are CM.
#  The rest is excluded for l >= 11)
