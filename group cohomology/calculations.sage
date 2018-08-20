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
