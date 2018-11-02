def solve_integer_problem_with_torsion(M, V, MT, VT, N, all=False):
    r"""
    Solves a linear problem in the integer where some relations are modulo N.

    INPUT:

    - ``M`` -- An m x n matrix with integer coefficients.
    - ``V`` -- A vector of length m with integer coefficients.
    - ``MT`` -- An mt x n matrix with coefficients in the integers
                or the integers modulo N.
    - ``VT`` -- A vector of length mt with coefficients in the
                integers or the integers modulo N.
    - ``N`` -- A non-negative integer giving the modulus of the problem.
    - ``all`` -- A boolean indicating whether a single solution or all
                 solutions to this problem should be returned.

    OUTPUT:
    
    A vector v of length n with integer coefficients such that
    M * v == V and MT * v == VT modulo N. If the parameter all is set
    to True, will return a tuple with as first entry such a vector v
    and as the second entry an n x n0 matrix A such that the map
    $y \mapsto v + A y$ maps n0 vectors with integer coefficients
    surjectively to the set of all solutions to the system.

    EXAMPLES::

    Let's find all pairs of integers that sum to 2018 and have the
    same last digit, i.e. congruent modulo 10::

        sage: M = matrix([1,1]); M
        [1 1]
        sage: V = vector([2018]); V
        (2018)
        sage: MT = matrix([1,-1]); MT
        [ 1 -1]
        sage: VT = vector([0]); VT
        (0)
        sage: N = 10
        sage: solve_integer_problem_with_torsion(M, V, MT, VT, N, all=True)
        (
                    [ 5]
        (-1, 2019), [-5]
        )
    """
    # Make sure the rings are right
    M = M.change_ring(ZZ)
    V = V.change_ring(ZZ)
    MT = MT.change_ring(Integers(N))
    VT = VT.change_ring(Integers(N))

    # Calculation + checking that a solution exists
    m, n = M.dimensions()
    mt, nt = MT.dimensions()
    M0 = block_matrix(2,2,[M, zero_matrix(m, mt), MT.change_ring(ZZ), diagonal_matrix([N]*mt)])
    V0 = vector(list(V) + list(VT.change_ring(ZZ)))
    MV0 = block_matrix(1,2,[M0,V0.column()])
    B = MV0.right_kernel().basis_matrix().transpose()
    g = B[-1][0]
    c = [1]
    for i in range(1, len(B[-1])):
        g, s, t = xgcd(g,B[-1][i])
        for j in range(len(c)):
            c[j] = c[j] * s
        c.append(t)
    if abs(g) != 1:
        raise ArithmeticError("The system does not have a solution.")
    v0 = -g * B * vector(c)
    v = v0[:n]
    if all:
        m0, n0 = MV0.dimensions()
        n0 = n0 - 1
        M0 = MV0.delete_columns([n0])
        A0 = M0.right_kernel().basis_matrix().transpose()
        A = A0[:n]
        return (v, A)
    else:
        return v
