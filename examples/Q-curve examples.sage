# Some examples:

### 2-isogeny parameter 1
# Singular

### 2-isogeny parameter 2
E = Qcurve_with_2_isogeny(2)
E = E.decomposable_twist()
K = E.decomposition_field()
E = E.twist(K(3))
f, twists = E.newform()
print f.Level() # 512 = 2^9
print f # q + a*q^3 - 2*a*q^5 - 4*q^7 - q^9 + a*q^11 + O(q^12)

### 2-isogeny parameter 3
E = Qcurve_with_2_isogeny(3)
E = E.decomposable_twist()
f, twists = E.newform()
print f.Level() # 1536 = 2^9 * 3
print f # q + (-a - 1)*q^3 + a*q^5 - 3*a*q^7 + (2*a - 1)*q^9 - 4*q^11 + O(q^12)

### 2-isogeny parameter 4
# Rational curve

### 2-isogeny parameter 5
E = Qcurve_with_2_isogeny(5)
E = E.decomposable_twist()
E = E.twist(QQ(3))
f, twists = E.newform()
print f.Level() # 800 = 2^5 * 5^2
print f # q + (-a + 1)*q^3 + (a + 1)*q^7 + a*q^9 + 4*a*q^11 + O(q^12)

### 2-isogeny parameter 6
E = Qcurve_with_2_isogeny(6)
E = E.decomposable_twist()
f, twists = E.newform()
print f.Level() # 7680 = 2^9 * 3 * 5
print f # q + 1/40*(-2*a^3 - 5*a^2 - 4*a - 30)*q^3 + 1/8*(a^2 + 6)*q^5 + 1/20*(a^3 + 22*a)*q^7 + 1/20*(-a^3 - 22*a + 20)*q^9 + 1/20*(a^3 + 2*a)*q^11 + O(q^12)

### 2-isogeny parameter -1
E = Qcurve_with_2_isogeny(-1)
E = E.decomposable_twist()
E = E.twist(QQ(3))
f, twists = E.newform()
print f.Level() # 512 = 2^9
print f # q + a*q^3 - 2*q^5 - 2*a*q^7 - q^9 - 3*a*q^11 + O(q^12)

### 2-isogeny parameter -2
E = Qcurve_with_2_isogeny(-2)
E = E.twist(QQ(3))
f, twists = E.newform()
print f.Level() # 1536 = 2^9 * 3
print f # q + q^3 + a*q^5 + a*q^7 + q^9 + 2*q^11 + O(q^12)

### 2-isogeny parameter -3
E = Qcurve_with_2_isogeny(-3)
E = E.decomposable_twist()
f, twists = E.newform()
print f.Level() # 1536 = 2^9 * 3 ? 96
print f # q + 1/16*(a^3 - 2*a^2 + 20*a - 24)*q^3 + 1/8*(a^3 + 28*a)*q^5 + 1/8*(a^3 + 20*a)*q^7 + 1/8*(-a^3 - 28*a + 8)*q^9 + 1/4*(a^2 + 12)*q^11 + O(q^12)

### 2-isogeny parameter -4
E = Qcurve_with_2_isogeny(-4)
E = E.decomposable_twist()
E = E.twist(QQ(3))
f, twists = E.newform()
print f.Level() # 1280 = 2^8 * 5
print f # q + a*q^3 + q^5 + 3*a*q^7 - q^9 + 4*a*q^11 + O(q^12)

### 2-isogeny parameter -5
E = Qcurve_with_2_isogeny(-5)
E = E.decomposable_twist()
E = E.twist(QQ(3))
f, twists = E.newform()
print f.Level() # 38400 = 2^9 * 3 * 5^2
print f # Memory problems

### 2-isogeny parameter -6
E = Qcurve_with_2_isogeny(-6)
E = E.decomposable_twist()
f, twists = E.newform()
print f.Level() # 10752 = 2^9 * 3 * 7
print f # q + 1/4*(-a^2 - 8*a - 14)*q^3 + 1/12*(a^3 + 12*a^2 + 58*a + 104)*q^5 + 1/24*(-a^3 - 12*a^2 - 58*a - 104)*q^7 + 1/2*(-a^2 - 8*a - 20)*q^9 - 4*q^11 + O(q^12)

### 3-isogeny parameter 0
# singular curve

### 3-isogeny parameter 1
# singular curve

### 3-isogeny parameter 2
# CM-curve

### 3-isogeny parameter 3
E = Qcurve_with_3_isogeny(3)
E = E.twist(QQ(2))
f, twists = E.newform()
print f.Level() # 3888 = 2^4 * 3^5
print f # q + a*q^5 + 6*q^11 + O(q^12)

### 3-isogeny parameter 4
# rational curve

### 3-isogeny parameter 5
E = Qcurve_with_3_isogeny(5)
E = E.decomposable_twist()
if (2^17).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(2))
f, twists = E.newform()
print f.Level() # 675 = 3^3 * 5^2
print f # q + 1/6*a^2*q^4 + a*q^7 + O(q^12)

### 3-isogeny parameter 6
E = Qcurve_with_3_isogeny(6)
E = E.decomposable_twist()
if (5^5).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(5))
f, twists = E.newform()
print f.Level() # 311040 = 2^8 * 3^5 * 5
print f # Memory problem

### 3-isogeny parameter -1
E = Qcurve_with_3_isogeny(-1)
E = E.decomposable_twist()
if (2^17).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(2))
f, twists = E.newform()
print f.Level() # 432 = 2^4 * 3^3
print f # q + 1/12*(a^3 + 30*a)*q^5 + 1/12*(a^3 + 18*a)*q^7 + 1/2*(-a^2 - 12)*q^11 + O(q^12)

### 3-isogeny parameter -2
E = Qcurve_with_3_isogeny(-2)
E = E.decomposable_twist()
f, twists = E.newform()
print f.Level() # 20736 = 2^8 * 3^4
print f # q + a*q^5 + 2*a*q^7 + O(q^12)

### 3-isogeny parameter -3
E = Qcurve_with_3_isogeny(-3)
E = E.decomposable_twist()
if (2^9).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(2))
f, twists = E.newform()
print f.Level() # 243 = 3^5
print f # q + a*q^2 + q^4 + 2*a*q^5 - q^7 - a*q^8 + 6*q^10 - 2*a*q^11 + O(q^12)

### 3-isogeny parameter -4
E = Qcurve_with_3_isogeny(-4)
E = E.decomposable_twist()
if (2^17).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(2))
if (5^5).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(5))
f, twists = E.newform()
print f.Level() # 2160 = 2^4 * 3^3 * 5
print f # q + 1/126*(-a^3 + 12*a^2 - 63*a + 124)*q^5 + 1/42*(-a^3 + 12*a^2 - 21*a - 44)*q^11 + O(q^12)

### 3-isogeny parameter -5
E = Qcurve_with_3_isogeny(-5)
E = E.decomposable_twist()
if (2^65).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(2))
f, twists = E.newform()
print f.Level() # 32400 = 2^4 * 3^4 * 5^2
print f # magma crash

### 3-isogeny parameter -6
E = Qcurve_with_3_isogeny(-6)
E = E.decomposable_twist()
if (7^5).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(5))
f, twists = E.newform()
print f.Level() # 435456 = 2^8 * 3^5 * 7
print f # Not tried, way too big

### 2,3-isogeny parameter 0
# Singular

### 2,3-isogeny parameter 1
# CM-curve

### 2,3-isogeny parameter 2
E = Qcurve_with_2_3_isogeny(2)
E = E.decomposable_twist() # Takes too long!
if (7^5).divides(E.conductor_restriction_of_scalars()):
    E = E.twist(QQ(5))
f, twists = E.newform()
print f.Level() # 
print f #

### 2,3-isogeny parameter 3
# CM-curve

### 2,3-isogeny parameter -1
# Not a Q-curve!

### 2,3-isogeny parameter -2
E = Qcurve_with_2_3_isogeny(-2)
E = E.complete_definition_twist([-2,-3])
E = E.decomposable_twist()
f, twists = E.newform()
print f.Level() # 2304 = 2^8 * 3^2
print f # q + 1/3335808*(-81*a^7 - 20*a^6 + 306*a^5 - 26616*a^4 - 152604*a^3 - 1406256*a^2 - 6279048*a + 2817760)*q^5 + 1/130560*(-4*a^7 + 3*a^6 - 1434*a^4 - 5904*a^3 - 70908*a^2 - 356800*a + 170568)*q^7 + 1/1308160*(5*a^7 - 24*a^6 - 246*a^5 + 2400*a^4 + 1244*a^3 - 58848*a^2 + 77816*a - 3780864)*q^11 + O(q^12)

### 2,3-isogeny parameter -3
