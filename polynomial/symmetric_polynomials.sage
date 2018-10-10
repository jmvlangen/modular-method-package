def polynomial_to_symmetric(polynomial, names=None):
    r"""
    Converts a symmetric polynomial to a polynomial in the
    the symmetric polynomials.

    INPUT:

    - ``polynomial`` -- A symmetric polynomial in any number
      of variables.
    - ``names`` -- A list of names for the symmetric
      polynomials. This list must have length n, where n
      is the number of variables in the parent of the given
      polynomial. By default will use the names s1, s2, s3,
      up to n.

    OUTPUT:
    
    A polynomial that when evaluated on the symmetric
    polynomials in which the given polynomial lives, returns
    the given polynomial.
    """
    R = polynomial.parent().base_ring()
    n = polynomial.parent().ngens()
    sym = SymmetricFunctions(R)
    f = sym.from_polynomial(polynomial)
    if names is None:
        names = ['s%s'%(i,) for i in range(1,n+1)]
    if len(names) != n:
        raise ValueError('The length of names is not %s'%(n,))
    S = PolynomialRing(R, names=names)
    gens = S.gens()
    g = sum(coeff * product((gens[i-1] if (i >= 1 and i <= n) else 0)
                            for i in partition)
            for partition, coeff in list(sym.elementary()(f)))
    return S(g)
