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

# Curve 1
# Y^2 = X^3 - 3*a*X - 2*c
a_invariants1 = [0, 0, 0, -3*a, -2*c]
E1 = FreyCurve(a_invariants1, condition=C)

# Curve2
# Y^2 = X^3 + 2*(sqrt(3) - 1)*(s - t)*X^2 + (2 - sqrt(3))*((s - t)^2 - 2*sqrt(3)*s*t)*X
K2.<sqrt3> = QuadraticField(3)
G2.<sigma> = K2.galois_group()
a_invariants2 = [0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0]
isogenies2 = {sigma^0: (QQ(1), 1), sigma^1: (-1 - sqrt3, 2)}
E2 = FreyQcurve(a_invariants2,
                isogenies=isogenies2,
                condition=C)

# The conductors
print ""
print E1.conductor()
# 64*3^n0*Rad_P( (-1728) * (s^4 - 4*s^3*t - 6*s^2*t^2 - 4*s*t^3 + t^4)^3 )
# where 
# n0 = 2 if ('s', 't') == (1, 2), (2, 1) mod 3
#      3 if ('s', 't') is 1 of 4 possibilities mod 3
# (agrees with article)
print ""
print E2.conductor()
# (64)*Rad_P( ((-960*sqrt3 + 1664)) * (s^2 + (2*sqrt3 - 2)*s*t + t^2) * (s^2 + (-2*sqrt3 - 2)*s*t + t^2)^2 )
# (agrees with article)

# Computing newforms with the multi-Frey method
print ""
nfs1 = apply_to_conditional_value(lambda x: list(x), E1.newform_candidates(algorithm='magma'))
nfs2 = apply_to_conditional_value(lambda x: list(x), E2.newform_candidates(algorithm='magma'))
nfs = eliminate_by_traces((E1, E2), (nfs1, nfs2), condition=C, primes=prime_range(5,20))
print nfs
# [(q + 4*q^5 + O(q^12), q + 1/2*(a + 2)*q^3 + a*q^5 + a*q^7 + (a - 1)*q^9 + 2*q^11 + O(q^12), 0), (q + 4*q^5 + O(q^12), q + 1/2*(-a + 2)*q^3 + a*q^5 + a*q^7 + (-a - 1)*q^9 + 2*q^11 + O(q^12), 0), (q + 4*q^5 + O(q^12), q + 1/2*(-a - 2)*q^3 + a*q^5 - a*q^7 + (a - 1)*q^9 - 2*q^11 + O(q^12), 0), (q + 4*q^5 + O(q^12), q + 1/2*(a - 2)*q^3 + a*q^5 - a*q^7 + (-a - 1)*q^9 - 2*q^11 + O(q^12), 0)]                   if ('s', 't') == (1, 2), (2, 1) mod 3
# [(q - 2*q^5 + 3*q^7 - 6*q^11 + O(q^12), q + 1/4*(a - 2)*q^3 + 1/2*(a - 8)*q^9 - 6*q^11 + O(q^12), 0), (q - 2*q^5 + 3*q^7 - 6*q^11 + O(q^12), q + 1/4*(-a + 2)*q^3 + 1/2*(a - 8)*q^9 + 6*q^11 + O(q^12), 0), (q - 2*q^5 - 3*q^7 + 6*q^11 + O(q^12), q + 1/4*(a - 2)*q^3 + 1/2*(a - 8)*q^9 - 6*q^11 + O(q^12), 0), (q - 2*q^5 - 3*q^7 + 6*q^11 + O(q^12), q + 1/4*(-a + 2)*q^3 + 1/2*(a - 8)*q^9 + 6*q^11 + O(q^12), 0)] if ('s', 't') is 1 of 4 possibilities mod 3
# (Agrees with article!)
