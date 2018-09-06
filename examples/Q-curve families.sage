# Curves with a 2-isogeny
def Qcurve_with_2_isogeny(t):
    K = field_with_root(QQ, t, names='sqrt_t')
    sqrt_t = sqrt(K(t))
    E = Qcurve([0,12,0,18*(sqrt_t + 1),0], guessed_degrees=[2])
    if E.has_cm():
        raise ValueError("Paramater %s gives a CM elliptic curve."%t)
    return E

# Curves with a 3-isogeny
def Qcurve_with_3_isogeny(t):
    K = field_with_root(QQ, t, names='sqrt_t')
    t = K(t)
    sqrt_t = sqrt(t)
    a4 = -3 * sqrt_t * (4 + 5*sqrt_t)
    a6 = 2 * sqrt_t * (2 + 14*sqrt_t + 11*t)
    E = Qcurve([a4,a6], guessed_degrees=[3])
    if E.has_cm():
        raise ValueError("Paramater %s gives a CM elliptic curve."%t)
    return E

# Curves with a 2 & 3-isogeny
def Qcurve_with_2_3_isogeny(t):
    # Parameter 1, 3, 9/2, 11/3, 25/6, 75/19, 289/72, 675/169, 1089/272, 31211/7803 are CM
    s = 1 + 2 * t
    K = composite_field(field_with_root(QQ, t, names='sqrt_t'), field_with_root(QQ, s, names='sqrt_s'))
    sqrt_t = sqrt(K(t))
    sqrt_s = sqrt(K(s))
    a4 = -6*s*t*(5 + 5*sqrt_s + 10*sqrt_t + 5*t + 2*sqrt_s*sqrt_t)
    a6 = 8 * (sqrt_s*sqrt_t)^3 * (1 + sqrt_t) * (7 + 15*sqrt_s + 14*sqrt_t + 7*t + 6*sqrt_s*sqrt_t)
    E = Qcurve([a4,a6], guessed_degrees=[2,3,6])
    if E.has_cm():
        raise ValueError("Paramater %s gives a CM elliptic curve."%t)
    return E
