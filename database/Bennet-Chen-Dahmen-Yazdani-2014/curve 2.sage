# Article: On the equation a^3 + b^(3n) = c^2
# Authors: M.A. Bennett, Imin Chen, Sander R. Dahmen, Soroosh Yazdani
# Journal: Acta Arith. 163 (2014), no. 4, 327--343
#
# Equation: a^3 + b^(3*l) = c^2
# Case: c odd
# Assumptions:
#   - a, b and c are coprime
#   - l > 7
#   - a = 2*(s^4 + 2*t*s^2 + 2*t^3*s + t^4)
#   - b^l = (s-t)^4 - 12*(s*t)^2
#   - c = 3*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)
#   - s != t (mod 2)
#   - s != t (mod 3)

# Variables
R.<s,t> = QQ[]
a = 2*(s^4 + 2*t*s^3 + 2*t^3*s + t^4)
bl = (s-t)^4 - 12*(s*t)^2
c = 3*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)

# Conditions
C = (CoprimeCondition(['s','t']) &
     ~CongruenceCondition(s-t, 2) &
     ~CongruenceCondition(s-t, 3) &
     PowerCondition(bl, 5))

# The curve
# Y^2 = X^3 + 2*(sqrt(3) - 1)*(s - t)*X^2 + (2 - sqrt(3))*((s - t)^2 - 2*sqrt(3)*s*t)*X
K.<sqrt3> = QuadraticField(3)
G.<sigma> = K.galois_group()
a_invariants = [0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0]
isogenies = {sigma^0: (QQ(1), 1), sigma^1: (-1 - sqrt3, 2)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               condition=C)

# The conductor
print ""
print E.conductor()
# (64)*Rad_P( ((-960*sqrt3 + 1664)) * (s^2 + (2*sqrt3 - 2)*s*t + t^2) * (s^2 + (-2*sqrt3 - 2)*s*t + t^2)^2 )
# (agrees with article)

# Newforms
print ""
# print E.newforms()
#
# (article: six forms remain of level 768)
