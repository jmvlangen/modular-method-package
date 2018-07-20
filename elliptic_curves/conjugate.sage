def conjugate_curve(E, sigma):
    """Conjugates the curve E by the automorphism sigma of a field."""
    return EllipticCurve([sigma(a) for a in E.a_invariants()])
