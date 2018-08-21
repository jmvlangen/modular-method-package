# Some examples:

### 2-isogeny parameter 2
E = Qcurve_with_2_isogeny(2)
E = E.decomposable_twist()
f, twists = E.newform()
print f # q + a*q^5 + 4*q^7 + 1/2*a*q^11 + O(q^12)

### 2-isogeny parameter 6
E = Qcurve_with_2_isogeny(6)
E = E.decomposable_twist()
f, twists = E.newform()
print f # q + 1/40*(-2*a^3 - 5*a^2 - 4*a - 30)*q^3 + 1/8*(a^2 + 6)*q^5 + 1/20*(a^3 + 22*a)*q^7 + 1/20*(-a^3 - 22*a + 20)*q^9 + 1/20*(a^3 + 2*a)*q^11 + O(q^12)

### 2-isogeny parameter -1
E = Qcurve_with_2_isogeny(-1)
E = E.decomposable_twist()
f, twists = E.newform()
print f # q - 2*q^5 + (-a - 2)*q^7 + 1/2*(-3*a - 6)*q^11 + O(q^12)

