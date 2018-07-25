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
        K = f.BaseField().sage()
        ap = f.Coefficient(p).sage()
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
            print factor, result
    if f.parent() == magma:
        return magma(result)
    else:
        return result

### Temporary stuff
candidates = range(8)
for p in prime_range(5, 100):
    i = 0
    print "-----------", p, "-----------"
    v = (f.Coefficient(p).sage()) * QQ(chi[1](p))
    print v
    for f1 in nfs0:
        ap = f1.Coefficient(p).sage()
        for s in ap.parent().galois_group():
            print s(ap), s(ap) == v
            if i in candidates and s(ap) != v:
                candidates.remove(i)
            i += 1
    print ""
print candidates
