# Article: Winding quotients and some variants of Fermat's Last Theorem
# Authors: Henri Darmon and Lo\"ic Merel
# Journal: J. reine angew. Math. 490 (1997), 81--100
# Source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
#
# Equation: a^l + b^l = c^3
# Case: a * b even
# Assumptions:
#   - a, b and c coprime
#   - l >= 7
#   - a odd (wlog)
#   - b even (wlog)

# Conditions
R.<bl, c> = ZZ[]
C = (CoprimeCondition([bl,c]) &
     PowerCondition(c^3 - bl, 7) & # c^3 - bl = a^l
     PowerCondition(bl, 7) & # bl = b^l
     CongruenceCondition(c^3 - bl -1, 2) & # wlog a^l odd
     CongruenceCondition(bl, 2)) # wlog

# The curve
# Y^2 + c*X*Y = X^3 - c^2*X^2 - 3/2*c*b^l*X + b^l*(a^l + 5/4*b^l))
a_invariants = [c, -c^2, 0, -3/2*c*bl, bl*((c^3 - bl) + 5/4*bl)]
E = FreyCurve(a_invariants, condition=C)

# Conductor
print ""
print E.conductor()
# 3^n0*Rad_P( (27) * bl * (c^3 - bl)^3 )
#  where 
# n0 = 2 if ('bl', 'c') is 1 of 710046 possibilities mod 2187
#      3 if ('bl', 'c') is 1 of 24 possibilities mod 9
#      1 if ('bl', 'c') is 1 of 1458 possibilities mod 2187
# (article: rad(a*b) if 3 divides a*b, 3^3*rad(a*b) otherwise)

# Newforms
print ""
print E.newforms()
# 
# (article: no solution if l == 1 (mod 3) and assuming Shimura-Tate)
