# Article: Winding quotients and some variants of Fermat's Last Theorem
# Authors: Henri Darmon and Lo\"ic Merel
# Journal: J. reine angew. Math. 490 (1997), 81--100
# Source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
#
# Equation: a^l + b^l = c^3
# Case: c even
# Assumptions:
#   - a, b and c coprime
#   - l >= 7
#   - c even

# Conditions
R.<bl, c> = ZZ[]
C = (CoprimeCondition([bl,c]) &
     PowerCondition(c^3 - bl, 7) & # c^3 - bl = a^l
     PowerCondition(bl, 7) & # bl = b^l
     CongruenceCondition(c, 2))

# The curve
# Y^2 + b^l*Y = X^3 - 3*((c/2)^3 + b^l)*(c/2)*X - (c/2)^3*(2*(c/2)^3 - 5*b^l)
a_invariants = [0, 0, bl, -3*((c/2)^3 + bl)*(c/2), -(c/2)^3*(2*(c/2)^3 - 5*bl)]
E = FreyCurve(a_invariants, condition=C)

# Conductor
print ""
N = E.conductor()
print N
# 3^n0*Rad_P( (27) * bl * (c^3 - bl)^3 )
#  where 
# n0 = 2 if ('bl', 'c') is 1 of 710046 possibilities mod 2187
#      3 if ('bl', 'c') is 1 of 24 possibilities mod 9
#      1 if ('bl', 'c') is 1 of 1458 possibilities mod 2187
# (article: rad(a*b) if 3 divides a*b, 3^3*rad(a*b) otherwise)

# Newforms
print ""
nfs = E.newforms(algorithm='magma')
print nfs
# 
# (article: no solution if l == 1 (mod 3) and assuming Shimura-Tate)
