# Article: On the equation a^3 + b^(3n) = c^2
# Authors: M.A. Bennett, Imin Chen, Sander R. Dahmen, Soroosh Yazdani
# Journal: Acta Arith. 163 (2014), no. 4, 327--343
#
# Equation: a^3 + b^(3*l) = c^2
# Case: c even
# Assumptions:
#   - a, b and c are coprime
#   - l > 7
#   - a = 3*s^4 + 6*t^2*s^2 - t^4
#   - b^l = -3*s^4 + 6*t^2*s^2 + t^4
#   - c = 6*s*t*(3*s^4 + t^4)
#   - t != 0 (mod 3)
#   - s != t (mod 2)

# Variables
R.<s,t> = QQ[]
a = 3*s^4 + 6*t^2*s^2 - t^4
bl = -3*s^4 + 6*t^2*s^2 + t^4
c = 6*s*t*(3*s^4 + t^4)

# Conditions
C = (CoprimeCondition(['s','t']) &
     ~CongruenceCondition(t, 3) &
     ~CongruenceCondition(s-t, 2) &
     PowerCondition(bl, 5))

# The curve
# Y^2 = X^3 - 3*a*X - 2*c
a_invariants = [0, 0, 0, -3*a, -2*c]
E = FreyQcurve(a_invariants, condition=C)

# The conductor
print ""
print E.conductor()
# 32*3^n0*Rad_P( (1728) * (3*s^4 - 6*s^2*t^2 - t^4)^3 )
#  where 
# n0 = 2 if ('s', 't') == (0, 1), (0, 2) mod 3
#      3 if ('s', 't') is 1 of 4 possibilities mod 3
# (agrees with article)

# Newforms
print ""
print E.newforms(algorithm='magma')
#
# (article: 1 form of level 288 and 2 forms of level 864 remain)
