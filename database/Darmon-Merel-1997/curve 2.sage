# Article: Winding quotients and some variants of Fermat's Last Theorem
# Authors: Henri Darmon and Lo\"ic Merel
# Journal: J. reine angew. Math. 490 (1997), 81--100
# Source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
#
# Equation: a^l + b^l = c^2
# Case: a*b even
# Assumptions:
#   - a, b and c coprime
#   - l >= 7
#   - a^l even (wlog)
#   - c == 1 (mod 4) (wlog)

# Conditions
R.<al, c> = QQ[]
C = (CoprimeCondition([al,c]) &
     PowerCondition(al, 7) & # al = a^l
     PowerCondition(c^2 - al, 7) & # c^2 - al = b^l
     CongruenceCondition(al, 2) & # wlog
     CongruenceCondition(c-1, 4)) # wlog

# The curve
# Y^2 + X*Y = X^3 + (c-1)/4 * X^2 + a^l / 2^6 X
a_invariants = [1, (c - 1)/4, 0, al / 2^6, 0]
E = FreyCurve(a_invariants, condition=C)

# Data
print ""
N = E.conductor()
print N
# Rad_P( (1/4096) * al^2 * (c^2 - al) )
# (agrees with article)

# Newforms
print ""
nfs = E.newforms(algorithm='magma')
print nfs
# 
# (article: none)
