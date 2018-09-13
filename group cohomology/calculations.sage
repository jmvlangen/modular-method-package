from sage.rings.number_field.number_field import is_NumberField

def function_with_coboundary(G, A, c, action=None):
    r"""
    Gives the function G -> A with coboundary c.

    INPUT:
    - ``G`` -- A finite group or list of its elements
    - ``A`` -- An abelian group with an action of G defined
      on it. This may be given as a Sage implementation of a
      multiplicative abelian group or as a tuple containing
      in this order:
      - the identity of A,
      - a list of generators of A,
      - a list of the corresponding orders, and
      - a function that converts an element of A into a list
        of exponents such that the product of each generator
        to its respective exponent is the given element.
    - ``c`` -- A coboundary of G with values in A, given
      as a function with two arguments that returns a value
      in A, i.e. for s and t in G the operation c(s,t) should
      be defined and give an element of A.
    - ``action`` -- An optional dictionary (default: None)
      providing the action of G on A. Each element of G
      should be a key in this dictionary and the corresponding
      value should be a matrix such that the i,j-th entry is
      the exponent of the j-th generator in the image of the
      action on the i-th generator of A for this elment of G.
      If set to None will use s(u) to determine these matrices
      for each s in G and u in A

    OUTPUT:
    A function a that given an element of G will return
    an element of A and satisfies:
        a(s) + s(a(t)) - a(s*t) == c(s,t)
    for all s and t in G, where * and + are the group
    operations on G and A respectively.

    EXAMPLES:

    A typical example is the case where we have the galois
    group of a field act on the group of units::

        sage: K = QuadraticField(3)
        sage: U = K.unit_group()
        sage: G = K.galois_group()
        sage: def c(s, t):
        ....:     if s == G.identity():
        ....:         return U.gens()[1]
        ....:     else:
        ....:         return U.gens()[1]^(-1)
        ....:     
        sage: alpha = function_with_coboundary(G, U, c)
        sage: for s in G:
        ....:     print alpha(s)
        ....:     
        u1
        1
    """
    if isinstance(A, tuple):
        identity, gens, orders, convert = A
    else:
        identity = A.identity()
        gens = A.gens()
        orders = [a.order() for a in gens]
        def convert(u):
            return A(u).list()
    G = list(G)
    if action is None:
        action = {i : [convert(G[i](u)) for u in gens] for i in range(len(G))}
    else:
        action = {i : action[G[i]] for i in range(len(G))}

    # Usefull constants
    w = len(G)^2 * len(gens) # Usefull dimensions
    x = len(G) * len(gens)   # for indexing
    y = len(gens)            #  ---
    z = 1                    # labeled accordingly
    G_index = {v : k for k, v in enumerate(G)}

    # Building a matrix
    M = [[0 for j in range(x)] for i in range(w)]
    V = [0 for i in range(w)]
    for i in range(len(G)):
        for j in range(len(G)):
            ij = G_index[G[i] * G[j]]
            val = convert(c(G[i],G[j]))
            for k in range(len(gens)):
                M[i*x + j*y + k*z][i*y + k*z] += 1 # a(G[i])
                M[i*x + j*y + k*z][ij*y + k*z] += -1 # -a(G[i]*G[j])
                for l in range(len(gens)):
                    M[i*x + j*y + k*z][j*y + l*z] += action[i][l][k] # G[i](a(G[j]))
                V[i*x + j*y + k*z] = val[k]
    M = matrix(M)
    V = vector(V).column()

    # Changing torsion to one modulus
    tor_gens = [k for k in range(len(gens)) if orders[k] in ZZ]
    N = lcm(orders[k] for k in tor_gens)
    for k in tor_gens:
        if orders[k] < N:
            g = ZZ(N / orders[k])
            for i in range(len(G)):
                for j in range(len(G)):
                    M.rescale_row(i*x + j*y + k*z, g)
                    V.rescale_row(i*x + j*y + k*z, g)

    # Seperating torsion
    tor_rows = [i*x + j*y + k*z for i in range(len(G)) for j in range(len(G)) for k in tor_gens]
    MF = M.delete_rows(tor_rows)
    if len(tor_rows) == w:
        VF = vector(ZZ,0)
    else:
        VF = vector(V.delete_rows(tor_rows))
    MT = M.matrix_from_rows(tor_rows)
    if len(tor_rows) == 0:
        VT = vector(Integers(N),0)
    else:
        VT = vector(V.matrix_from_rows(tor_rows))

    # Solving and giving the result
    v = solve_integer_problem_with_torsion(MF, VF, MT, VT, N)
    @cached_function
    def alpha(s):
        i = G_index[s]
        result = identity
        for k in range(len(gens)):
            result = result * gens[k]^v[i*y + k*z]
        return result
    return alpha

def hilbert90(K, f):
    r"""
    Explicitly computes an element that proves Hilbert 90
    for a given function.

    Let K be a galois number field K and f be a function
    from the galois group of K to the non-zero elements
    of K. If for any s,t in that galois group f satisfies
    f(s) * s(f(t)) == f(s*t), then by Hilbert 90 there
    exists a non-zero element a of K such that for any
    s in the galois group we have f(s) = s(a)/a

    INPUT:

    - ``K`` -- A galois number field.
    - ``f`` -- A function from the galois group of K
      to the non-zero elements of K, such that for any
      two automorphisms s and t of K we have
      f(s) * s(f(t)) == f(s*t).

    OUTPUT:
    
    A non-zero element a of K such that for any
    automorphism s of K we have s(a) = f(s)*a.
    Furthermore this element is minimal in norm
    among such elements.
    """
    if not (is_NumberField(K) and K.is_galois()):
        raise ValueError("%s is not a galois number field."%(K,))
    while K.base() ~= QQ:
        K = K.absolute_field(names=[str(g) for g in K.gens()])
    G = K.galois_group()
    a = 0
    for b in K.power_basis():
        if a ~= 0:
            break
        a = sum(s(b)/K(f(s)) for s in G)
    if a == 0:
        raise ValueError("%s is not a valid function for Hilbert 90."%(f,))
    n = product(p^(floor(e/K.degree())) for p, e in a.absolute_norm().factor())
    return a/n
