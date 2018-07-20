def function_with_coboundary(G, A, c):
    r"""
    Gives the function G -> A with coboundary c.

    INPUT:
    - ``G`` -- A finite group
    - ``A`` -- An abelian group with an action of G defined
               on it, i.e. for each s in G and a in A the
               operation s(a) should be defined and give an
               element of A.
    - ``c`` -- A coboundary of G with values in A, given
               as a function with two arguments that returns
               a value in A, i.e. for s and t in G the
               operation c(s,t) should be defined and give
               an element of A.

    OUTPUT:
    A function a that given an element of G will return
    an element of A and satisfies:
        a(s) + s(a(t)) - a(s*t) == c(s,t)
    for all s and t in G, where * and + are the group
    operations on G and A respectively.
    """
    gens = A.gens()
    G = list(G)
    G_on_gens = {i : [A(G[i](u)).list() for u in gens] for i in range(len(G))}

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
            val = A(c(G[i],G[j])).list()
            for k in range(len(gens)):
                M[i*x + j*y + k*z][i*y + k*z] += 1 # a(G[i])
                M[i*x + j*y + k*z][ij*y + k*z] += -1 # -a(G[i]*G[j])
                for l in range(len(gens)):
                    M[i*x + j*y + k*z][j*y + l*z] += G_on_gens[i][l][k] # G[i](a(G[j]))
                V[i*x + j*y + k*z] = val[k]
    M = matrix(M)
    V = vector(V).column()

    # Changing torsion to one modulus
    tor_gens = [k for k in range(len(gens)) if gens[k].order() in ZZ]
    N = lcm(gens[k].order() for k in tor_gens)
    for k in tor_gens:
        if gens[k].order() < N:
            g = ZZ(N / gens[k].order())
            for i in range(len(G)):
                for j in range(len(G)):
                    M.rescale_row(i*x + j*y + k*z, g)
                    V.rescale_row(i*x + j*y + k*z, g)

    # Seperating torsion
    tor_rows = [i*x + j*y + k*z for i in range(len(G)) for j in range(len(G)) for k in tor_gens]
    MF = M.delete_rows(tor_rows)
    VF = vector(V.delete_rows(tor_rows))
    MT = M.matrix_from_rows(tor_rows)
    VT = vector(V.matrix_from_rows(tor_rows))

    # Solving and giving the result
    v = solve_integer_problem_with_torsion(MF, VF, MT, VT, N)
    @cached_function
    def alpha(s):
        i = G_index[s]
        result = A.identity()
        for k in range(len(gens)):
            result = result * gens[k]^v[i*y + k*z]
        return result
    return alpha
