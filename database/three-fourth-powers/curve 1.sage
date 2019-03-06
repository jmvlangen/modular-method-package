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
# Y^2 = X^3 - 40*b*X^2 - 20*(Sqrt(30)*a^2 - 10*b^2 + 2*Sqrt(30)*y^2)*X
K.<w> = QuadraticField(30)
G.<sigma> = K.galois_group()
a_invariants = [0, -40*b, 0, -20*(w*a^2 - 10*b^2 + 2*w*b^2),0]
L.<sqrtm2> = QuadraticField(-2)
isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm2, 2)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               condition=C)

# Twisting
E = E.decomposable_twist()

# Conductor
N = E.conductor()
print ""
print N
#

# The levels of the newforms
levels = E.newform_levels()
print ""
print levels
# [(15360, 15360, 76800, 76800), (76800, 76800, 15360, 15360)]

# Newforms
nfs = E.newforms(algorithm='file', path='tmp/15360_full_magma.nfs')
print ""
print len(nfs)
# 
