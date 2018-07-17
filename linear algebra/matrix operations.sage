def smallest_entry_in_column(M, j, N=0, i0=0, give_gcd=False):
    r"""
    Finds the row with the smallest entry in a given column.

    This method works on integer matrices or matrices defined
    over the integers modulo some number. The smallest value
    is determined by the gcd of the number and the modulus.
    Note that for modulus 0 this corresponds to working over
    $\Z$ and using the standard absolute value.

    INPUT:

    - ``M`` -- A matrix with integer entries or entries in
               the ring of integers modulo a given N
    - ``j`` -- The column in which one wants to find the
               smallest element. Must be an integer between
               0 and M.dimensions()[1]
    - ``N`` -- A non-negative integer (default=0) giving the
               modulus with which to work.
    - ``i0`` -- A non-negative integer (default=0) indicating
                the least row to be checked. All rows before
                the i0-th row will be ignored. Setting this value
                to at least M.dimensions()[0] will return 0 by
                default.
    - ``give_gcd`` -- A boolean (default=False) indicating
                      whether the corresponding gcd should also
                      be returned.

    OUTPUT:
    
    Returns the index of the row of M of which the j-th entry has
    the smallest possible gcd with N possible. Note that gcd 0 will
    be considered to be the largest possible gcd. If no entry was
    found, e.g. if i0 >= M.dimensions()[0], then i0 will be returned
    as the index. If the flag give_gcd is set to True, the returned
    value will be a tuple containing the index as the first entry
    and the corresponding gcd as its second value. Note that this
    gcd will be N if i0 >= M.dimensions()[0].

    EXAMPLES:
    
    For modulus 0 finds the smallest value ::

        sage: M = matrix([[1,2,0],[0,3,0]]); M
        [1 2 0]
        [0 3 0]
        sage: smallest_entry_in_column(M, 1)
        0

    whilst for modulus 4 it finds the smallest gcd ::
    
        sage: M = matrix([[1,2,0],[0,3,0]]); M
        [1 2 0]
        [0 3 0]
        sage: smallest_entry_in_column(M, 2, N=4, give_gcd=True)
        (1, 1)
    
    Note that it considers 0 as the biggest value ::
    
        sage: M = matrix([[1,2,0],[0,3,0]]); M
        [1 2 0]
        [0 3 0]
        sage: smallest_entry_in_column(M, 0, give_gcd=True)
        (0, 1)

    Forcing it to look only below a certain row we can force
    if to choose 0 as the smallest ::
    
        sage: M = matrix([[1,2,0],[0,3,0]]); M
        [1 2 0]
        [0 3 0]
        sage: smallest_entry_in_column(M, 0, i0=1, give_gcd=True)
        (1, 0)

    Choosing a too big i0 will always result in i0 and 0 ::

        sage: M = matrix([[1,2,0],[0,3,0]]); M
        [1 2 0]
        [0 3 0]
        sage: smallest_entry_in_column(M, 2, i0=5, give_gcd=True)
        (5, 0)

    """
    gcd_min = N
    i_min = i0
    for i in range(i0, M.dimensions()[0]):
        g = gcd(ZZ(M[i][j]), N)
        if g != N and (gcd_min == N or g < gcd_min):
            gcd_min = g
            i_min = i
    if give_gcd:
        return i_min, gcd_min
    else:
        return i_min

def lower_zero_rows(M):
    r"""
    Swaps rows in a matrix to move zero rows to the bottom.

    Swaps rows in a given matrix M until all zero rows are
    at the bottom. Does not necessarily preserve the order
    of the non-zero rows.

    INPUT:
    
    - ``M`` -- A matrix

    OUPUT:

    None

    EXAMPLE::

        sage: M = matrix([[1,0,0],[0,0,0],[0,2,3]]); M
        [1 0 0]
        [0 0 0]
        [0 2 3]
        sage: lower_zero_rows(M); M
        [1 0 0]
        [0 2 3]
        [0 0 0]

    """
    m,n = M.dimensions()
    for i0 in range(m):
        if M[i0] == 0: # Zero row
            flag = True
            for i1 in range(i0+1, m):
                if M[i1] != 0: # first non-zero row thereafter
                    flag = False
                    break
            if flag: # done
                break
            M.swap_rows(i0, i1)

def eliminate_zero_rows(M):
    r"""
    Removes all non-zero rows from a matrix.

    INPUT:

    - ``M`` -- A matrix

    OUTPUT:
    
    The matrix M without all the zero rows.

    EXAMPLE::

        sage: M = matrix([[1,0,0],[0,0,0],[0,2,3]]); M
        [1 0 0]
        [0 0 0]
        [0 2 3]
        sage: eliminate_zero_rows(M)
        [1 0 0]
        [0 2 3]

    """
    m, n = M.dimensions()
    ls = [i for i in range(m) if M[i] == 0]
    return M.delete_rows(ls)

def unit_for_scaling(a, b, N=0):
    r"""
    Finds a unit u such that a * u = b.

    The unit is found modulo N. If no such unit exists,
    will give an arithmetic error.

    INPUT:
    
    - ``a`` - An integer or an integer modulo N
    - ``b`` - An integer or an integer modulo N such
              that gcd(a,N) == gcd(b,N)
    - ``N`` - A non-negative integer (default=0) that
              indicates the modulus to work with

    OUTPUT:

    An integer u that is a unit modulo N such
    that a * u == b modulo N.

    EXAMPLES:

    Modulo 0 it simply gives the right sign ::

        sage: unit_for_scaling(5,-5)
        -1
        sage: unit_for_scaling(3,3)
        1

    Modulo we can test some of the results ::

        sage: u = unit_for_scaling(294, 6, N=1728); u
        241
        sage: u * 294 % 1728
        6
        sage: gcd(u, 1728)
        1

    Note that a and b must differ by a unit modulo N ::

        sage: unit_for_scaling(294, 1248, N = 1728)
        ArithmeticError: There is not a unit u modulo 1728 such that 294 * u == 1248 modulo 1728
    """
    g = gcd(ZZ(a),N)
    if g != gcd(ZZ(b),N):
        raise ArithmeticError("There is not a unit u modulo %s such that %s * u == %s modulo %s"%(N, a, b, N))
    if N == 0:
        return sign(ZZ(a) * ZZ(b))
    M = N / g
    aM = Integers(M)(ZZ(a)/g)
    bM = Integers(M)(ZZ(b)/g)
    uM = bM / aM
    for u in range(ZZ(uM), N, M):
        if gcd(u, N) == 1:
            break
    return u

def normalize_row(M, i, j, val=0, N=0):
    r"""
    Rescales a row based on a given entry in that row.

    For a given matrix M, multiplies row i by a unit modulo N
    such that a given entry j in that row takes on the
    value val, which by default will be set to the gcd of
    that entry with N.

    INPUT:
    
    - ``M`` -- A matrix with integer coefficients or
               coefficients in the integers modulo N.
    - ``i`` -- A non-negative integer smaller than
               M.dimensions()[0], indicating the row
               to normalize.
    - ``j`` -- A non-negative integer smaller than
               M.dimensions()[1], indicating the entry
               to normalize.
    - ``val`` -- The value to normalize to, which by
                 default (when set to 0) will be the
                 gcd of M[i][j] and N.
    - ``N`` -- A non-negative integer (default=0)
               indicating the modulus to work with.

    OUPUT:

    None

    EXAMPLES:

    Modulo 0 we can only rescale by plus or minus 1 ::
    
        sage: M = matrix([[1,-2,0],[0,3,0]]); M
        [ 1 -2  0]
        [ 0  3  0]
        sage: normalize_row(M, 0, 1); M
        [-1  2  0]
        [ 0  3  0]
        sage: normalize_row(M, 0, 0); M
        [ 1 -2  0]
        [ 0  3  0]
        sage: normalize_row(M, 1, 1, val=1); M
        ArithmeticError: There is not a unit u modulo 0 such that 3 * u == 1 modulo 0

    Modulo some number there is more possibilities ::
    
        sage: M = matrix([[20,21,0],[0,22,0]]).change_ring(Integers(24)); M
        [20 21  0]
        [ 0 22  0]
        sage: normalize_row(M, 0, 0, N=24); M
        [ 4  9  0]
        [ 0 22  0]
        sage: normalize_row(M, 0, 1, N=24); M
        [20  3  0]
        [ 0 22  0]
        sage: normalize_row(M, 1, 1, val=10, N=24); M
        [20  3  0]
        [ 0 10  0]
    """
    if val == 0:
        val = gcd(ZZ(M[i][j]), N)
    M.rescale_row(i, unit_for_scaling(M[i][j], val, N))
    
def eliminate_entry_with_row(M, i, j):
    r"""
    Subtracts a row from other rows to kill their j-th entry.

    Given a matrix M and a row i of that matrix, will attempt
    to substract row i from the other rows in M such that
    j-th coefficient in these rows becomes smaller than the
    j-th coefficient in row i.

    INPUT:
    
    - ``M`` -- A matrix with coefficients in the integers or
               a quotient ring thereof.
    - ``i`` -- A non-negative integer smaller than
               M.dimensions()[0] indicating which row to
               use to eliminate entries in the others.
    - ``j`` -- A non-negative integer smalle than
               M.dimensions()[1] indicating the entry to
               be eliminated.

    OUTPUT:

    None

    EXAMPLES:

    Can kill both positive and negative values ::

        sage: M = matrix([[2,3],[1,-1],[-3,4]]); M
        [ 2  3]
        [ 1 -1]
        [-3  4]
        sage: eliminate_entry_with_row(M, 1, 0); M
        [ 0  5]
        [ 1 -1]
        [ 0  1]

    Note that it will not eliminate smaller values ::
    
        sage: M = matrix([[2,3],[1,-1],[-3,4]]); M
        [ 2  3]
        [ 1 -1]
        [-3  4]
        sage: eliminate_entry_with_row(M, 0, 0); M
        [ 2  3]
        [ 1 -1]
        [ 1 10]

    """
    if M[i][j] == 0:
        raise ArithmeticError("Can not eliminate rows using 0.")
    for i0 in range(M.dimensions()[0]):
        if i0 != i:
            c = -floor(ZZ(M[i0][j]) / ZZ(M[i][j]))
            M.add_multiple_of_row(i0, i, c)

def minimal_echelon_form(M, N=0):
    r"""
    Puts a matrix in echelon form without the zero rows.

    INPUT:

    - ``M`` -- A matrix with coefficients in the integers
               or in the integers modulo N.
    - ``N`` -- A non-negative integer (default=0)
               indicating the modulus we should work with.

    OUTPUT:
    
    A matrix which is an echelon form of M without all
    the zero rows. Here an echelon form means that for
    every pivot:
    - it is normalized, i.e. its gcd with N is equal
      to itself
    - all entries below it are zero
    - all entries above it are smaller
    
    NOTE:

    If N=0 and we work over the integers the echolon
    form might contain some rows of which the entries
    have gcd > 1. The user should keep this into
    account if using this method for calculating a
    kernel.

    EXAMPLE::

        sage: M = matrix([[-4, 1, -3],[2, -4, -5]]); M
        [-4  1 -3]
        [ 2 -4 -5]
        sage: minimal_echelon_form(M)
        [ 2  3  8]
        [ 0  7 13]
        sage: minimal_echelon_form(M.change_ring(Integers(24)), N=24)
        [ 2  0 23]
        [ 0  1 19]
    """
    i0 = 0
    for j in range(M.dimensions()[1]):
        i, g = smallest_entry_in_column(M, j, N=N, i0=i0, give_gcd=True)
        if g != N:
            g0 = g + 1 # Make it loop the first time!
            while g < g0: # Loop to make sure we definitely have the smallest g
                M.swap_rows(i0, i)
                normalize_row(M, i0, j, val=g, N=N)
                eliminate_entry_with_row(M, i0, j)
                M = eliminate_zero_rows(M)
                g0 = g
                i, g = smallest_entry_in_column(M, j, N=N, i0=i0, give_gcd=True)
            i0 += 1
    return M

def right_kernel(M, N=0):
    r"""
    Gives the right kernel of a matrix M as a matrix.

    Gives the biggest matrix A with linearly independent
    columns such that M * A = 0

    INPUT:
    
    -- ``M`` - A matrix with coefficients in the integers
               or the integers modulo N
    -- ``N`` - A non-negative integer indicating the
               modulus to work with.

    OUTPUT:
    
    A matrix A with linearly independent columns over the
    integers modulo N that satisfies M * A == 0 modulo N.
    Any other such matrix has at most as many columns
    as A.
    """
    M0 = minimal_echelon_form(M, N=N)
    m, n = M0.dimensions()
    if N == 0: # Maybe once isn't enough?
        for i in range(m):
            g = gcd(list(M0[i]))
            M0.rescale_row(i, 1/g)
        M0 = minimal_echelon_form(M)
    pivot_ls = []
    for i in range(m):
        for j in range(n):
            if M0[i][j] != 0:
                ls.append(j)
                break
    A = []
    for j in range(n):
        if j not in pivot_ls:
            v = [0]*n
            v[j] = -1
            ## TODO
