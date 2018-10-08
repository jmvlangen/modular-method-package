# Article: Winding quotients and some variants of Fermat's Last Theorem
# Authors: Henri Darmon and Lo\"ic Merel
# Journal: J. reine angew. Math. 490 (1997), 81--100
# Source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
#
# Equation: a^l + b^l = c^2
# Case: a*b odd
# Assumptions:
#   - a, b and c coprime
#   - l >= 7
#   - a^l == -1 (mod 4) (wlog)

# Conditions
R.<al, c> = QQ[]
C = (CoprimeCondition([al,c]) &
     PowerCondition(al, 7) & # al = a^l
     PowerCondition(c^2 - al, 7) & # c^2 - al = b^l
     CongruenceCondition(al + 1, 4)) # wlog

# The curve
# Y^2 = X^3 + 2 * c * X^2 + a^l X
a_invariants = [0, 2*c, 0, al, 0]
E = FreyCurve(a_invariants, condition=C)

# Conductor
print ""
print E.conductor()
# 32*Rad_P( 2^6 * ap^2 * (c^2 - ap) )
# (agrees with article)

# Newforms
print ""
print E.newforms()
# 
# (article: eliminated by other means (CM))
