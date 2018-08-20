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
