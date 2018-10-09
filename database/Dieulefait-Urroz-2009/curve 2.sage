# Article: Solving Fermat-type equations via modular Q-curves over polyquadratic fields
# Authors: Luis Dieulefait and Jorge Jim\'enez Urroz
# Journal: J. Reine Angew. Math. 633 (2009), 183–195
# Source (arxiv): https://arxiv.org/abs/math/0611663
#
# Equation: a^4 + 3*b^2 = c^l
# Assumptions:
#   - a, b and c are coprime
#   - l > 131

# Variables
R.<a,b> = QQ[]
d = 3
cl = a^4 + d*b^2

# Conditions
C = (CoprimeCondition([a,b]) &
     PowerCondition(cl, 5))

# The curve
# Y^2 = X^3 + 4*a*X^2 + 2*(a^2 + sqrt(-d)*b)*X
K.<r> = QuadraticField(-d)
L.<sqrtm2> = QuadraticField(-2)
G.<sigma> = K.galois_group()
a_invariants = [0, 4*a, 0, 2*(a^2 + r*b), 0]
isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm2, 2)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               condition=C)

# Conductor
print ""
print E.conductor()
# (16)*Rad_P( ((-233610467125568/15*rsqrtm2rzeta0^7 + 583050224176000*rsqrtm2rzeta0^6 - 32212260985250176/15*rsqrtm2rzeta0^5 + 10494904035168000*rsqrtm2rzeta0^4 + 114728482001153792/15*rsqrtm2rzeta0^3 - 34983013450560000*rsqrtm2rzeta0^2 - 583609468822355456/15*rsqrtm2rzeta0 + 222330073930210304)) * (a^2 + (1/4800*rsqrtm2rzeta0^7 + 7/800*rsqrtm2rzeta0^5 + 63/400*rsqrtm2rzeta0^3 + 169/600*rsqrtm2rzeta0)*b) * (a^2 + (-1/4800*rsqrtm2rzeta0^7 - 7/800*rsqrtm2rzeta0^5 - 63/400*rsqrtm2rzeta0^3 - 169/600*rsqrtm2rzeta0)*b)^2 )
# (agees with article)

# Twisting
E = E.decomposable_twist()
# This is not the same twist as in the article!

# Newform
print ""
# print E.newforms()
#
# (article: no newforms agree for l > 131)
