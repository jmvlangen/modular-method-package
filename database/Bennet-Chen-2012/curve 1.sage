# Article: Multi-Frey Q-curves and the Diophantine equation a^2 + b^6 = c^n
# Authors: Michael A. Bennett and Imin Chen
# Journal: Algebra & Number Theorey, volume 6 (2012), no. 4
# Source: http://people.math.sfu.ca/~ichen/pub/BeCh2.pdf
#
# Equation: a^2 + b^2 = c^l
# Assumptions:
#   - a, b and c are coprime

# Variables
R.<a,b> = QQ[]
cl = a^2 + b^6

# Conditions
C = (CoprimeCondition([a,b]) &
     PowerCondition(cl, 5))

# The curve
K.<i> = QuadraticField(-1)
L.<sqrtm3> = QuadraticField(-3)
G.<sigma> = K.galois_group()
a_invariants = [0, 0, 0, -3*(5*b^3 + 4*a*i)*b, 2*(11*b^6 + 14*i*b^3*a - 2*a^2)]
isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm3 , 3)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               condition=C)

# Twisting
E = E.decomposable_twist()
# Agrees with the twist found in the article

# Conductor
print ""
N = E.conductor()
print N
# (4)*(1/4*izeta0^2 + 1)^n0*Rad_P( ((18195840*izeta0^3 - 145566720*izeta0 + 252129024)) * (b^3 + (-1/8*izeta0^3)*a) * (b^3 + (1/8*izeta0^3)*a)^3 )
#  where 
# n0 = 0 if ('a', 'b') is 1 of 24 possibilities mod 9
#      4 if ('a', 'b') is 1 of 48 possibilities mod 9
# (agrees with article)

# Newforms
print ""
nfs = E.newforms(algorithm='magma')
print nfs
#
# (article: the only form at level 48 is CM and 2 forms at level 432 are CM.
#  The rest is excluded for l >= 11)
