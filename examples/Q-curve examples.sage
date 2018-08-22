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
print f # q + (a + 1)*q^3 + a*q^5 + 3*a*q^7 + (2*a - 1)*q^9 + 4*q^11 + O(q^12)

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
print f.Level() # 1536 = 2^9 * 3
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

