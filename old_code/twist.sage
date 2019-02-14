def compute_possible_twists(K, P):
    r"""
    Computes twists that could change the conductor of an elliptic curve.

    Given an elliptic curve over a number field, the conductor exponent
    at a given prime might change if the curve is twisted. This function
    computes a list of elements of the number field, such that all these
    changes in the conductor at a given prime can be achieved by twisting
    by one of the elements in this list.

    INPUT:
    
    - ``K`` -- A number field, possibly $\\Q$, over which the elliptic
               curve would be defined.
    - ``P`` -- A prime ideal of K.

    OUTPUT:
    A list of elements of K. Every change in the conductor exponent
    at P of an elliptic curve defined over K, corresponds to at least
    one twist of the curve by an element in this list. The trivial
    twist (1) is always part of this list.

    EXAMPLES:

    Well known results for $\\Q$::

        sage: compute_possible_twists(QQ, 2)
        [1, 3, 2, 6]
        sage: compute_possible_twists(QQ, 5)
        [1, 5]

    Works over bigger number fields::

        sage: K = CyclotomicField(7)
        sage: compute_possible_twists(K, K.prime_above(2))
        [1,
        -zeta7^5 + zeta7^4 - zeta7^3,
        -zeta7^5 + zeta7^4 + zeta7^3,
        -zeta7^5 - zeta7^4 + zeta7^3,
        zeta7^5 - zeta7^4 + zeta7^3,
        zeta7^5 + zeta7^4 - zeta7^3,
        -zeta7^5 - zeta7^4 - zeta7^3,
        zeta7^5 + zeta7^4 + zeta7^3,
        zeta7^3 + zeta7 + 1,
        2*zeta7^5 + 2*zeta7^4 + zeta7^3 + 2*zeta7^2 + zeta7 + 3,
        2*zeta7^4 + zeta7^3 - zeta7 + 1,
        -2*zeta7^5 + zeta7^3 - zeta7 - 1,
        -2*zeta7^5 - 2*zeta7^4 - zeta7^3 - 2*zeta7^2 - zeta7 - 3,
        2*zeta7^5 - zeta7^3 + zeta7 + 1,
        zeta7^3 + 2*zeta7^2 + zeta7 + 1,
        -zeta7^3 - 2*zeta7^2 - zeta7 - 1]
        sage: compute_possible_twists(K, K.prime_above(3))
        [1, 3]

    """
    base = pAdicBase(K, P)
    pi = base.uniformizer()
    if base.characteristic() != 2:
        return [K(1), pi]
    else:
        result =  _compute_possible_unit_twists(base)
        return result + [c * pi for c in result]

def _compute_possible_unit_twists(base):
    r"""
    Returns the different twists of an elliptic curve by units.

    INPUT:
    - ``base`` -- A pAdicBase object that contains the field and prime

    OUTPUT:
    Similar to the output of :func: compute_possible_twists,
    but only computes those that are units modulo the prime.
    """
    K = base.number_field()
    T_unram = _compute_unramified_unit_tree(base)
    T_ram = T_unram.complement()
    T_ram.root().children.get((K(0),)).remove()
    k = T_unram.root().minimum_full_level()
    unram = [node.quotient_tuple()[0] for node in T_unram.nodes_at_level(k)]
    ram = [node.quotient_tuple()[0] for node in T_ram.nodes_at_level(k)]
    result = [K(1)]
    while len(ram) > 0:
        val = ram[0]
        result.append(val.lift())
        for val2 in unram:
            if val * val2 in ram:
                ram.remove(val * val2)
            else:
                "Warning: possibly a mistake in `_compute_possible_unit_twists`"
    return result

def _tree_of_unit_squares_mod(base, mod):
    r"""
    Compute the squares modulo some power of a prime in a local ring.

    INPUT:
    - ``base`` -- a pAdicBase object that describes the local ring.
    - ``mod`` -- the power of the prime modulo which we want to
                 consider squares.

    OUTPUT:
    A pAdicTree that contains all numbers in the local ring,
    that are squares modulo the prime to the power mod.
    """
    K = base.number_field()
    k = base.valuation(2)
    n = min(mod - k, floor((mod + 1)/2))
    T = pAdicTree(variables='x', pAdics=base)
    T.root().children.get((K(0),)).remove()
    T2 = pAdicTree(variables='x', pAdics=base, full=False)
    for node in T.nodes_at_level(n):
        val = node.representative()[0]^2
        node = T2.root()
        for coeff in base.power_series(val, mod):
            coeffs = (coeff,)
            if not node.children.contains(coeffs):
                node.children.add(pAdicNode(parent=node,
                                            coefficients=coeffs))
            node = node.children.get(coeffs)
        if not node.is_full():
            node.children = pAdicNodeCollection_inverted(node)
    return T2
        
def _compute_unramified_unit_tree(base):
    r"""
    Computes the tree of units of a local ring that give unramified
    quadratic extensions.

    INPUT:
    - ``base`` -- A pAdicBase object that describes the local ring

    OUTPUT:
    A pAdicTree that contains all units of the local ring, such that
    the extension by the square root of these remains unramified.
    """
    k = base.valuation(2)
    T = _tree_of_unit_squares_mod(base, 2*k)
    for i in reversed(range(0,2*k,2)):
        for node in T.nodes_at_level(i):
            for coeffs in base.representatives():
                if not node.children.contains(coeffs) \
                and not (i == 0 and coeffs[0] == 0):
                    node.children.add(pAdicNode(parent=node,
                                                coefficients=coeffs,
                                                full=True))
    return T
