### Init
R.<a,b> = QQ[]
f = (a-b)^4 + a^4 + (a+b)^4

### Phase 1: Limiting units
S.<x> = QQ[]
L.<v> = NumberField(f(x,1))
# Over this field we have f = 3*(a + v*b)*(a - v*b)*(a^2 + (v^2 + 4)*b^2)
h = a + v*b
# Outside primes above 2, 3 and 5, h and its galois conjugates
# in the field over which f splits are coprime.
# Using the knowledge that c is not divisible by 2, 3 and 5
# and the fact that v has valuation -1 at the prime above 3 only
# we see that the valuation of (a + v*b) must be
#  0 at primes above 2 and 5
#  -1 at the prime above 3
# Since P3^4 = 3, where P3 is the prime above 3, we find that
#   h = 3 * u * ch^5 (*)
# for u a unit and ch in K
# 
# Since the valuation of h at a prime is only negative at the prime above 3,
# the valuation of ch is non-negative at all primes except the prime above 3,
# Furthermore at the prime above 3 the valuation of ch is precisely -1
# This implies that 3*ch is in the ring of integers, hence
R4.<s1, s2, s3, s4> = QQ[]
ch = 1/3 * sum(product(term) for term in zip(R4.gens(), L.integral_basis()))
UL = L.unit_group()
# Note that UL has two generators of which UL.0 has order 2 and is thus an 5-th power
# Therefore we only have to consider as units UL.1^i for i = 0, ..., 5
# This gives us the cases for the RHS in (*)
RHS = [L(3) * UL.1^i * ch^5 for i in range(5)]
# Given a basis of L containing 1 and v, we can write equation (*) as the a system
# of equations over the rationals
B = [1, v, v^2, v^3]
RHSQ = [polynomial_split_on_basis(hi, B) for hi in RHS]
# For each case i the system of equations is
#   a = RHSLQ[i][0]
#   b = RHSLQ[i][1]
#   0 = RHSLQ[i][2]
#   0 = RHSLQ[i][3]
# Note that the equations are integral outside 3,
# hence we can consider them modulo primes not equal to 3
RHS5 = [[tmp2.change_ring(GF(5)) for tmp2 in tmp1] for tmp1 in RHSQ]
# For each case we have a way of writing the pair a,b
# using the fact that the last two entries of RHS5 should be zero
print (RHS5[0][0],               RHS5[0][1] -  RHS5[0][3]) # (a,b) = (-,0) (mod 5)
print (RHS5[1][0] -  RHS5[1][2], RHS5[1][1] +  RHS5[1][3]) # (a,b) = (0,0) (mod 5)
print (RHS5[2][0] +  RHS5[2][2], RHS5[2][1])               # (a,b) = (0,0) (mod 5)
print (RHS5[3][0],               RHS5[3][1] +2*RHS5[3][3]) # (a,b) = (0,0) (mod 5)
print (RHS5[4][0] +2*RHS5[4][2], RHS5[4][1] +  RHS5[4][3]) # (a,b) = (0,-) (mod 5)
# Since a and b are coprime, only the first and last case remain

### Phase 2: Exluding the remaining cases through hyperelliptic curves
# L has a natural automorphism sending v to -v
s = L.hom([-v], L)
# Now g = h * s(h) has coefficients in the field K = Q(-v^2)
# which is a quadratic field that also contains the other factor of f
K.<w>, KtoL = L.subfield(-v^2)
g = a + w*b
# From (*) it follows that g satisfies
#   g = 9 * u' * cg^5 (**)
# where cg = ch * s(ch) is an element of K
# and u' = u * s(u) is a unit in K
# Note that f = g * t(g), where t is the non-trivial automorphism on K
# Again we know that cg is integral except possible at the prime above 3
# Using that f = g * t(g), equation (**) and the valuation of w,
# we find that the valuation of cg at the prime above 3 is -1,
# hence 3*cg is integral in K.
R2.<t1,t2> = QQ[]
cg = 1/3*sum(product(term) for term in zip(R2.gens(), K.integral_basis()))
# For the units we have two possibilities coming from those obtained for K
U = [K(u*s(u)) for u in [UL.1^i for i in [0,4]]]
# This gives us the possibilities of the RHS of (**)
RHS = [K(9) * u * cg^5 for u in U]
# Note that we can change the unit by a fifth power of a unit in K
# without changing all the properties we wrote above of cg.
# We will do so to get more suitable hyperelliptic curves later on
UK = K.unit_group()
RHS[0] = RHS[0] * UK.1^5
# Again we can use this to rewrite (**) to equations over Q
B = [1, w]
RHSQ = [polynomial_split_on_basis(gi, B) for gi in RHS]
# The equations are for case i
#   a^2 = RHSQ[i][0]
#   b^2 = RHSQ[i][1]
# Looking at these equations we can see that t1 must be divisible
# by 3 as a^2 and b^2 are integral, hence we can substitute 3*t1 for t1
RHSQ = [[poly(3*t1, t2) for poly in case] for case in RHSQ]
# We can exclude the case t2 = 0 as for that case a^2 and b^2 are not coprime
print [[poly(t1, 0).change_ring(GF(3)) for poly in case] for case in RHSQ]
# Multiplying the two equations for a case and dividing by t2^10 gives
# us the affine part of a hyperelliptic curve in terms of x = (t1/t2) and y = (a*b/t2)
Hs = [HyperellipticCurve(product(case)(x, 1)) for case in RHSQ]
Hm = [magma(H) for H in Hs]
# Note that both curves have odd degree models as they contain a rational point
HmOdd = [H.HasOddDegreeModel(nvals=2)[1] for H in Hm]
# The rational points correspond to the linear factors of RHSQ[0][1] and RHSQ[1][0] respectively
print RHSQ[0][1].factor() # (5) * (6*t1 + t2) * (419049*t1^4 + 295272*t1^3*t2 + 78036*t1^2*t2^2 + 9168*t1*t2^3 + 404*t2^4)
print RHSQ[1][0].factor() # (-5) * (23*t1 + 88*t2) * (201580749*t1^4 + 3084714072*t1^3*t2 + 17701580436*t1^2*t2^2 + 45146766768*t1*t2^3 + 43178995204*t2^4)
# Since t1 and t2 must be coprime, as a^2 and b^2 are
# these correspond to the points
#   (t1, t2) = (1, -6), (-1, 6) for the first case
#   (t1, t2) = (88, -23), (-88, 23) for the second case
# which give rise to the points (a^2, b^2):
print (RHSQ[0][0](1, -6), RHSQ[0][1](1, -6)) # (a^2, b^2) = (9, 0)
print (RHSQ[0][0](-1, 6), RHSQ[0][1](-1, 6)) # (a^2, b^2) = (-9, 0)
print (RHSQ[1][0](88, -23), RHSQ[1][1](88, -23)) # (a^2, b^2) = (0, 4)
print (RHSQ[1][0](-88, 23), RHSQ[1][1](-88, 23)) # (a^2, b^2) = (0, -4)
# Hence they do not correspond to solutions
# We can find bounds on the elements of their jacobians
Jm = [H.Jacobian() for H in Hm]
JmOdd = [H.Jacobian() for H in HmOdd]
# Note that we use an odd degree model for the first case for the rankbound on
# its jacobian as it gives a lower bound.
print [JmOdd[0].RankBound(), Jm[1].RankBound()] # [0, 0]
print [J.TorsionBound(50) for J in Jm] # [4, 4]
# Hence both jacobians contain at most 4 elements, hence only torsion elements
# The torsion maps injectively into the reduction modulo a good prime.
print [H.BadPrimes() for H in Hm] # [[ 2, 3, 5 ], [ 2, 3, 5 ]]
# So we can glance the torsion structure from looking modulo 7
Jm7 = [H.ChangeRing(GF(7)).Jacobian() for H in Hm]
print [ZZ(len(J.RationalPoints())).factor() for J in Jm7] # [2^2 * 5 * 11^2, 2^2 * 5 * 11^2]
# So the torsion must be Z/2Z, Z/2Z x Z/2Z or Z/4Z in both cases
print [(lcm([g.Order().sage() for g in J.AbelianGroup().Generators()])).factor() for J in Jm7] # (2 * 5 * 11^2, 2 * 5 * 11^2)
# So all torsion points on the Jacobians must be two torsion points
# These we can compute explicitly, since both curves have an odd degree model
print [J.TwoTorsionSubgroup().Order() for J in JmOdd] # [2, 2]
# So both Jacobians contain only two points
# One of each arises from the rational point on the corresponding hyperelliptic curve
# that does not correspond to a solution
# The other point corresponds to (the points on) the divisor
#   RHSQ[0][0](x, 1) = y = 0 for the first case
#   RHSQ[1][1](x, 1) = y = 0 for the second case
# Hence there exist no points on the hyperelliptic curves corresponding to solutions
# Hence there are no solutions!
