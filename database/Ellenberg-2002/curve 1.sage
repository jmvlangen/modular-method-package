# Article: Galois representations attached to Q-curves and
#          the generalized Fermat equation A^4 + B^2 = C^p
# Authors: Jordan S. Ellenberg
# Journal: Amer. J. Math. 126(4), 763--787 (2004)
# Source: http://www.math.wisc.edu/~ellenber/A4B2Cp.pdf
#
# Equation: a^4 + b^4 = c^l
# Assumptions:
#   - a, b and c coprime
#   - l >= 211
#   - b != 1 (mod 4) (wlog)

# Conditions
R.<a,b> = QQ[]
C = (CoprimeCondition([a,b]) &
     PowerCondition(a^4 + b^2, 3) & # a^4 + b^2 = c^p
     ~CongruenceCondition(b-1, 4)) # b != 1 mod 4

# The curve
# Y^2 = X^3 + 2*(1 + i)*a*X^2 + (b + i*a^2)*X
K.<i> = QuadraticField(-1)
G.<s> = K.galois_group()
a_invariants = [0, 2*(1+i) * a, 0 , b + i*a^2, 0]
isogenies={s^0: (QQ(1), 1), s^1: (1+i, 2)}
E = FreyQcurve(a_invariants,
               isogenies=isogenies,
               parameter_ring=ZZ,
               condition=C)

# Conductor
print ""
N = E.conductor()
print N
# Fractional ideal (i + 1)^n0*Rad_P( ((-64*i)) * (a^2 + (i)*b) * (a^2 + (-i)*b)^2 )
#  where 
# n0 = 12 if ('a', 'b') == (1, 0) mod 2
#      6  if ('a', 'b') == (0, 3), (2, 3) mod 4
# (agrees with article)

# Newforms
print ""
nfs = E.newforms(algorithm='magma')
print nfs
#
# (article: None if l >= 211)
