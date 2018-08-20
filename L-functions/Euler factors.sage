def Euler_factor_modular_form(f, p, twists=[1]):
    r"""
    Computes the Euler factor of a modular form.

    INPUT:

    - ``f`` -- A modular form. This may be a magma
      or a sage object.
    - ``p`` -- A prime number.

    OUTPUT:

    The p-th Euler factor of the L-Series, either
    as a magma polynomial if f was a magma object
    or as a sage polynomial otherwise.
    """
    if f.parent() == magma:
        ap = f.Coefficient(p).sage()
        K = ap.parent()
        epsp = f.DirichletCharacter()(p).sage()
    else:
        K = f.base_ring()
        ap = f.coefficient(p)
        epsp = f.character()(p)
    R.<T> = K[]
    result = R.one()
    if twists[0] not in ZZ:
        twists = [K(chi(p)) for chi in twists]
    for chip in twists:
        for F in K.galois_group():
            factor = (1 - F(chip * ap) * T + F(chip^2 * epsp * p) * T^2)
            result = result * factor
    if f.parent() == magma:
        return magma(result)
    else:
        return result

def _Euler_factor_elliptic_curve_local(E, P, R):
    r"""
    Computes the Euler factor of an elliptic curve corresponding to a given prime.

    INPUT:

    - ``E`` -- An elliptic curve over a number field.
    - ``P`` -- A prime in the number field over which E is defined.
    - ``R`` -- A polynomial ring in one variable over QQ

    OUTPUT:

    A polynomial in R that is the factor of the Euler factor of E that
    corresponds to the prime P.
    """
    NP = P.norm()
    p, e = NP.factor()[0]
    x = R.gen()^e
    if E.has_additive_reduction(P):
        return R.one()
    elif E.has_nonsplit_multiplicative_reduction(P):
        return R.one() + x
    elif E.has_split_multiplicative_reduction(P):
        return R.one() - x
    else:
        aP = NP + 1 - E.reduction(P).count_points()
        return R.one() - aP * x + NP * x^2

def _Euler_factor_elliptic_curve_local_QQ(E, p, R):
    r"""
    Computes the Euler factor of an elliptic curve over $\Q$.

    INPUT:

    - ``E`` -- An elliptic curve over $\Q$
    - ``p`` -- A prime number.
    - ``R`` -- A polynomial ring in one variable over QQ

    OUTPUT:

    A polynomial in R that is the Euler factor of E that
    corresponds to the prime number p.
    """
    x = R.gen()
    if E.has_additive_reduction(p):
        return R.one()
    elif E.has_nonsplit_multiplicative_reduction(p):
        return R.one() + x
    elif E.has_split_multiplicative_reduction(p):
        return R.one() - x
    else:
        ap = p + 1 - E.reduction(p).count_points()
        return R.one() - ap * x + p * x^2

def Euler_factor_elliptic_curve(E, p):
    r"""
    Compute the p-th Euler factor of an elliptic curve.

    INPUT:
    
    - ``E`` -- An elliptic curve defined over a number field.
    - ``p`` -- A prime number.

    OUTPUT:

    The Euler factor corresponding to p of the L-series
    of the given elliptic curve.

    """
    R.<T> = QQ[]
    K = E.base_ring()
    if K == QQ:
        return _Euler_factor_elliptic_curve_local(E, p, R)
    else:
        primes = K.primes_above(p)
        return product(_Euler_factor_elliptic_curve_local(E, P, R) for P in primes)
