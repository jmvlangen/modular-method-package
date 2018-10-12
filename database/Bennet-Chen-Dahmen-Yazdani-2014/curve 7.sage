# Article: On the equation a^3 + b^(3n) = c^2
# Authors: M.A. Bennett, Imin Chen, Sander R. Dahmen, Soroosh Yazdani
# Journal: Acta Arith. 163 (2014), no. 4, 327--343
#
# Equation: a^3 + b^(3*l) = c^2
# Case: c even
# Assumptions:
#   - a, b and c are coprime
#   - l > 7
#   - a = -3*s^4 + 6*t^2*s^2 + t^4
#   - b^l = (t^2 + 3*s^2)^2 - 12*s^4
#   - c = 6*s*t*(3*s^4 + t^4)
#   - s and t coprime
#   - t != 0 (mod 3)
#   - s != t (mod 2)

# Variables
R.<s,t> = QQ[]
a = 3*s^4 + 6*t^2*s^2 - t^4
bl = (t^2 + 3*s^2)^2 - 12*s^4
c = 6*s*t*(3*s^4 + t^4)

# Conditions
C = (CoprimeCondition([s,t]) &
     ~CongruenceCondition(t, 3) &
     ~CongruenceCondition(s-t, 2) &
     PowerCondition(bl, 5))

# The curve
# Y^2 = X^3 + 12*(sqrt(3) - 1)*s*X^2 - 3*sqrt(3)*(sqrt(3) - 1)^2*(t^2 + (2*sqrt(3)+3)*s^2)*X
K.<sqrt3> = QuadraticField(3)
G.<sigma> = K.galois_group()
a_invariants = [0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3+3)*s^2), 0]
isogenies={sigma^0: (QQ(1), 1), sigma^1: (-1 - sqrt3, 2)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               condition=C)

# The conductor
print ""
print E.conductor()
# (192)*Rad_P( ((-1492992*sqrt3 + 2612736)) * (s^2 + (-2/3*sqrt3 - 1)*t^2) * (s^2 + (2/3*sqrt3 - 1)*t^2)^2 )
# (agrees with article)

# Newforms
print ""
print E.newforms(algorithm='magma')
#
# (article: 4 forms remain of level 2304)
