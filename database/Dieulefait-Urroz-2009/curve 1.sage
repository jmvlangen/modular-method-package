# Article: Solving Fermat-type equations via modular Q-curves over polyquadratic fields
# Authors: Luis Dieulefait and Jorge Jim\'enez Urroz
# Journal: J. Reine Angew. Math. 633 (2009), 183–195
# Source (arxiv): https://arxiv.org/abs/math/0611663
#
# Equation: a^4 + 2*b^2 = c^l
# Assumptions:
#   - a, b and c are coprime
#   - l > 349

# Variables
R.<a,b> = QQ[]
d = 2
cl = a^4 + d*b^2

# Conditions
C = (CoprimeCondition([a,b]) &
     PowerCondition(cl, 5))

# The curve
# Y^2 = X^3 + 4*a*X^2 + 2*(a^2 + sqrt(-d)*b)*X
K.<r> = QuadraticField(-d)
G.<sigma> = K.galois_group()
a_invariants = [0, 4*a, 0, 2*(a^2 + r*b), 0]
isogenies = {sigma^0: (QQ(1), 1), sigma^1: (r, 2)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               condition=C)

# Conductor
print ""
print E.conductor()
# (r)^n0*Rad_P( (512) * (a^2 + (-r)*b) * (a^2 + (r)*b)^2 )
#  where 
# n0 = 12 if ('a', 'b') == (1, 1) mod 2
#      10 if ('a', 'b') == (1, 0) mod 2
# (agees with article)

# Newform
print ""
# print E.newforms()
#
# (article: no newforms agree for l > 349)
