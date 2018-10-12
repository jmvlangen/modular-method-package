# Article: Winding quotients and some variants of Fermat's Last Theorem
# Authors: Henri Darmon and Lo\"ic Merel
# Journal: J. reine angew. Math. 490 (1997), 81--100
# Source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
#
# Equation: a^l + b^l = 2*c^l
# Assumptions:
#   - a, b and c coprime
#   - l >= 7
#   - a^l == -1 (mod 4) (wlog)

# Conditions
R.<al, cl> = QQ[]
C = (CoprimeCondition([al,cl]) &
     PowerCondition(al, 7) & # al = a^l
     PowerCondition(2*cl - al, 7) & # 2*cl - al = b^l
     PowerCondition(cl, 7) & # cl = c^l
     CongruenceCondition(al + 1, 4)) # wlog

# The curve
# Y^2 = X*(X - a^l)*(X - 2*c^l)
a_invariants = [0, (-al)+(-2*cl), 0, (-al)*(-2*cl), 0]
E = FreyCurve(a_invariants, condition=C)

# Data
print ""
print E.conductor()
# 2^n0*Rad_P( 2^6 * cl^2 * al^2 * (al - 2*cl)^2 )
#  where 
# n0 = 5 if ('al', 'cl') == (3, 1), (3, 3) mod 4
#      1 if ('al', 'cl') is 1 of 8 possibilities mod 32
# (agrees with article)

# Newforms
print ""
print E.newforms(algorithm='magma')
# 
# (article: none if 2 divides a*b*c)
