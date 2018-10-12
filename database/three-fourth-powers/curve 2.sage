# Equation: (a - b)^4 + a^4 + (a + b)^4 = c^l
# Assumptions:
#   - a, b and c are coprime
#   - l >= 5

# Variables
R.<a,b> = QQ[]
cl = (a - b)^4 + a^4 + (a + b)^4

# Conditions
C = (CoprimeCondition([a,b]) &
     PowerCondition(cl, 5)) # cl = c^l

# The curve
# Y^2 = X^3 - 60*a*X^2  - 30*(-15*a^2 + 3*Sqrt(30)*a^2 + Sqrt(30)*b^2)*X
K.<w> = QuadraticField(30)
G.<sigma> = K.galois_group()
a_invariants = [0, -60*a, 0, -30*(-15*a^2 + 3*w*a^2 + w*b^2),0]
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
#

# The levels of the newforms
print ""
print E._newform_levels()
# [(23040, 23040, 115200, 115200), (115200, 115200, 23040, 23040)] if ('a', 'b') == (1, 0) mod 2 and ('a', 'b') == (1, 0) mod 2
# []                                                               if ('a', 'b') == (1, 0) mod 2 and ('a', 'b') == (1, 1) mod 2
# []                                                               if ('a', 'b') == (1, 1) mod 2 and ('a', 'b') == (1, 0) mod 2
# [(11520, 11520, 57600, 57600), (57600, 57600, 11520, 11520)]     if ('a', 'b') == (1, 1) mod 2 and ('a', 'b') == (1, 1) mod 2

# Newforms
print ""
print E.newforms(algorithm='magma')
# 
