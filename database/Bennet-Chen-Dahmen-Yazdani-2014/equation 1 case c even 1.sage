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

# Curve 1
# Y^2 = X^3 - 3*a*X - 2*c
a_invariants1 = [0, 0, 0, -3*a, -2*c]
E1 = FreyCurve(a_invariants1, condition=C)

# Auxiliary field
K.<sqrt3> = QuadraticField(3)
G.<sigma> = K.galois_group()

# Curve 2
# Y^2 = X^3 + 4*(sqrt(3) - 1)*t*X^2 - (sqrt(3) - 1)^2*(sqrt(3)*s^2 + (-2 - sqrt(3))*t^2)*X
a_invariants2 = [0, 4*(sqrt3 - 1)*t, 0, -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 - sqrt3)*t^2), 0]
isogenies2={sigma^0: (QQ(1), 1), sigma^1: (-1 - sqrt3, 2)}
E2 = FreyQcurve(a_invariants2,
                isogenies=isogenies2,
                condition=C)

# Curve 3
# Y^2 = X^3 + 12*(sqrt(3) - 1)*s*X^2 - 3*sqrt(3)*(sqrt(3) - 1)^2*(t^2 + (2*sqrt(3)+3)*s^2)*X
a_invariants3 = [0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3+3)*s^2), 0]
isogenies3={sigma^0: (QQ(1), 1), sigma^1: (-1 - sqrt3, 2)}
E3 = FreyQcurve(a_invariants3,
               isogenies=isogenies3,
               condition=C)

# The conductors
print ""
print E1.conductor()
# 32*3^n0*Rad_P( (1728) * (3*s^4 - 6*s^2*t^2 - t^4)^3 )
#  where 
# n0 = 2 if ('s', 't') == (0, 1), (0, 2) mod 3
#      3 if ('s', 't') is 1 of 4 possibilities mod 3
# (agrees with article)
print ""
print E2.conductor()
# (64)*Rad_P( ((39936*sqrt3 - 69120)) * (s^2 + (2/3*sqrt3 - 1)*t^2) * (s^2 + (-2/3*sqrt3 - 1)*t^2)^2 )
# (agrees with article)
print ""
print E3.conductor()
# (192)*Rad_P( ((-1492992*sqrt3 + 2612736)) * (s^2 + (-2/3*sqrt3 - 1)*t^2) * (s^2 + (2/3*sqrt3 - 1)*t^2)^2 )
# (agrees with article)

# Computing newforms with the multi-Frey method
print ""
nfs1 = apply_to_conditional_value(lambda x: list(x), E1.newform_candidates(algorithm='magma'))
nfs2 = apply_to_conditional_value(lambda x: list(x), E2.newform_candidates(algorithm='magma'))
nfs3 = apply_to_conditional_value(lambda x: list(x), E3.newform_candidates(algorithm='magma'))
nfs = eliminate_newforms_by_trace((E1, E2, E3), (nfs1, nfs2, nfs3), condition=C, primes=prime_range(5,20))
print nfs
# 
# (Agrees with article!)
