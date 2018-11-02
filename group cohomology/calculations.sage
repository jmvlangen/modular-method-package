from sage.rings.number_field.number_field import is_NumberField

def function_with_coboundary(G, A, c, action=None):
    r"""
    Gives the function G -> A with coboundary c.

    INPUT:
    - ``G`` -- A finite group 
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
    if action is None:
        action = {s : [convert(s(u)) for u in gens] for s in G}
    action = {s : matrix(action[s]).transpose() for s in action}

    # Constructing the important relations
    relations = [(s, t) for s in G for t in G]
    G_without_1 = [s for s in G if s != G.identity()]
    # Removing some relations
    relations.remove((G.identity(), G.identity()))
    for s in G_without_1:
        if (G.identity(), s) in relations:
            relations.remove((G.identity(), s)) # Linearly dependent with (1,1)
        if (s, G.identity()) in relations:
            relations.remove((s, G.identity())) # Linearly dependent with (1,1)
        if (s^2, s) in relations and (s, s^2) in relations:
            relations.remove((s^2, s)) # Linearly dependent with (s, s^2)
        for t in G_without_1:
            for w in G_without_1:
                if (not (s == t and t == w) and
                    (t, w) in relations and (s*t, w) in relations and
                    (s, t*w) in relations and (s, t) in relations):
                    relations.remove((t,w)) # These four are distinct and linearly dependent

    # The variables
    m = len(G.gens()) # Number of generators of G
    n = len(gens) # Number of generators of A
    alpha0 = {G.gens()[i]: block_matrix(1, m+1,
                                        [(identity_matrix(n) if j == i else zero_matrix(n))
                                         for j in range(m)] +
                                        [vector([0]*n).column()])
              for i in range(m)}
    def c_to_matrix(s,t):
        return block_matrix([[zero_matrix(n) for j in range(m)] +
                             [vector(convert(c(s, t))).column()]])
    # We always have alpha(1) = c(1, 1)
    alpha0[G.identity()] = c_to_matrix(G.identity(), G.identity())
    while len(alpha0) < len(G):
        keys = alpha0.keys()
        for s in keys:
            for t in keys:
                if s*t not in alpha0:
                    alpha0[s*t] = alpha0[s] + (action[s] * alpha0[t]) - c_to_matrix(s, t)
                    if (s,t) in relations:
                        relations.remove((s,t))

    # Building a matrix:
    MV = block_matrix(len(relations), 1, [alpha0[s] + action[s]*alpha0[t] - alpha0[s*t] - c_to_matrix(s, t)
                                         for (s, t) in relations])

    # Changing torsion to one modulus
    tor_gens = [k for k in range(n) if orders[k] in ZZ]
    N = lcm(orders[k] for k in tor_gens)
    for k in tor_gens:
        if orders[k] < N:
            g = ZZ(N / orders[k])
            for i in range(len(relations)):
                MV.rescale_row(i*n + k, g)

    # Seperating torsion
    tor_rows = [i*n + k for i in range(len(relations)) for k in tor_gens]
    MVF = MV.delete_rows(tor_rows)
    MVT = MV.matrix_from_rows(tor_rows)

    # Extracting the distinct matrices
    MF = MVF[:,:-1]
    VF = -MVF[:,-1]
    if VF.dimensions()[0] == 0:
        VF = vector([])
    else:
        VF = vector(VF)
    MT = MVT[:,:-1]
    VT = -MVT[:,-1]
    if VT.dimensions()[0] == 0:
        VT = vector([])
    else:
        VT = vector(VT)

    # Solving and giving the result
    v = solve_integer_problem_with_torsion(MF, VF, MT, VT, N)
    @cached_function
    def alpha(s):
        result = identity
        val = alpha0[s] * block_matrix(2,1,[v.column(), identity_matrix(1)])
        for k in range(n):
            result = result * gens[k]^val[k][0]
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
    Furthermore this element has the minimal
    integral norm among such elements.
    """
    if not (is_NumberField(K) and K.is_galois()):
        raise ValueError("%s is not a galois number field."%(K,))
    while K.base() != QQ:
        K = K.absolute_field(names=[str(g) for g in K.gens()])
    G = K.galois_group()
    if not all(f(s) * s(f(t)) == f(s*t) for s in G for t in G):
        raise ValueError("%s is not a valid function for Hilbert 90."%(f,))
    a = 0
    for b in K.power_basis():
        if a != 0:
            break
        a = sum(s(b)/K(f(s)) for s in G)
    if a == 0:
        raise ArithmeticError("Disproved Hilbert 90")
    I = a.numerator_ideal()
    J = a.denominator_ideal()
    d = I.smallest_integer()
    m = (K.ideal(d)*J*I^(-1)).smallest_integer()    
    n = d / m # Biggest rational such that a is an integral
              # element times that element.
    return a/n
