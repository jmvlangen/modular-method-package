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
# Y^2 = X^3 - 3*a*X - 2*c
a_invariants = [0, 0, 0, -3*a, -2*c]
E = FreyCurve(a_invariants, condition=C)

# The conductor
print ""
print E.conductor()
# 64*3^n0*Rad_P( (-1728) * (s^4 - 4*s^3*t - 6*s^2*t^2 - 4*s*t^3 + t^4)^3 )
# where 
# n0 = 2 if ('s', 't') == (1, 2), (2, 1) mod 3
#      3 if ('s', 't') is 1 of 4 possibilities mod 3
# (agrees with article)

# Newforms
print ""
print E.newforms()
#
# (article: 1 form of level 576 and 2 forms of level 1728 remain)
