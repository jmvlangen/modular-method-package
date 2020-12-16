r""" Some functions to compute with p-adic extension fields

AUTHORS:

- Joey van Langen (2020-12-12): initial version

"""
from modular_method.padics.pAdic_base import pAdicBase
from modular_method.padics.pAdic_tree import pAdicTree
from modular_method.padics.pAdic_solver import find_pAdic_roots

def _quot_comp(ls, subls, I):
    r"""Compute classes of `subls` in `ls` modulo `I`"""
    result = []
    done = []
    for a1 in ls:
        if a1 not in done:
            extra = [a3 for a3 in ls if any(a3 - a2*a1 in I for a2 in subls)]
            done.extend(extra)
            result.append(extra)
    return result

def ramified_quadratic_extensions(pAdics=None, ring=None, prime=None,
                                  inert=False, unramified=False,
                                  trivial=False):
    r"""Compute representatives for classes of ramified quadratic
    extensions.

    Given a local field $K$ each quadratic extension is of the form
    $K(\sqrt{\gamma}$ for some $\gamma \in K^*$. There is a subgroup
    $U \subseteq K^*$ consisting of those $\gamma \in K^*$ for which
    $K(\sqrt{\gamma}) / K$ is unramified. This $U$ in particular
    contains $(K^*)^2$ which corresponds to the trivial extensions. The
    group $K^*/U$ is therefore an $\ZZ/2\ZZ$-vector space.

    Note $K^*$ is isomorphic to $\pi^\ZZ \ZZ_K^*$ for a uniformizer
    $\pi$ of the maximal ideal of the ring of integers $\ZZ_K$ of
    $K$. By Hensel lifting every element of $\ZZ_K$ is a square if and
    only if it is a square modulo a high enough power of
    $\pi$. Therefore $K^* / (K^*)^2$ and also $K^* / U$ are finite and
    we can find representatives thereof in $\ZZ_K$ moduo a high enough
    power of $\pi$.

    INPUT:

    - ``pAdics`` -- A pAdicBase object (default: pAdicBase(ring,
      prime)) that defines the local field $K$

    - ``ring`` -- A ring (default: prime.ring()) that can be used to
      define a pAdicBase object. Will be ignored if `pAdics` is
      defined.

    - ``prime`` -- A prime ideal (default: pAdics.prime_ideal()) of a
      number field that defines the local field. Will be ignored if
      pAdics is defined.

    - ``inert`` -- A boolean value (default: False). If set to true
      will return additional information about extensions in which the
      prime is inert.

    - ``unramified`` -- A boolean value (default: False). If set to true
      will return additional information about unramified extensions

    - ``trivial`` -- A boolean value (default: False). If set to true
      will return additional information about trivial extensions

    OUTPUT:

    A list of elements of the number field corresponding to the given
    p-adics, one for each class in $K^*/U$.

    If `inert` was set to True will return an additional list that is
    a list of representatives from the number field corresponding to
    the given p-adics, one for each class of $U/(K^*)^2$.

    If `unramified` was set to True will return an additional list of
    elements of the number field corresponding to the given
    p-adics. The set $V$ of all elements in $\ZZ_K^*$ that reduce to
    one of the elements in this list modulo $4$ satisfies $U = \pi^{2
    \ZZ} V$.

    If `trivial` was set to True will return an additional list of
    elements of the number field corresponding to the given
    p-adics. The set $W$ of all elements in $\ZZ_K^*$ that reduce to
    one of the elements in this list modulo $4 \pi$ satisfies $(K^*)^2
    = \pi^{2 \ZZ} W$.

    EXAMPLES:

    Some simple examples over localizations of the rationals::

        sage: from modular_method.padics.pAdic_extensions import ramified_quadratic_extensions
        sage: ramified_quadratic_extensions(ring=QQ, prime=2)
        [1, 3, 2, 6]
        sage: ramified_quadratic_extensions(ring=QQ, prime=3)
        [1, 3]
        sage: ramified_quadratic_extensions(ring=QQ, prime=5)
        [1, 5]

    It also works for localizations of number fields::

        sage: from modular_method.padics.pAdic_extensions import ramified_quadratic_extensions
        sage: K = QuadraticField(3)
        sage: ramified_quadratic_extensions(prime=K.prime_above(2))
        [1, 6*a + 11, a + 2, 7*a + 12, a + 1, 17*a + 29, 3*a + 5, 19*a + 33]
        sage: ramified_quadratic_extensions(prime=K.prime_above(3))
        [1, a]
        sage: ramified_quadratic_extensions(prime=K.prime_above(5))
        [1, 5]

    The different quantities you can compute correspond to one another::

        sage: from modular_method.padics.pAdic_extensions import ramified_quadratic_extensions
        sage: K = QuadraticField(-2)
        sage: P2 = K.prime_above(2)
        sage: ram, inert, unram, triv = ramified_quadratic_extensions(prime=K.prime_above(2),
        ....: inert=True, unramified=True, trivial=True)
        sage: len(ram) * len(inert) * len(triv) == (P2^5).norm()
        True
        sage: len(ram) * len(unram) == (P2^4).norm()
        True
        sage: all(any(a - b in P2^4 for b in unram) for a in triv)
        True
        sage: all(any(a - b in P2^4 for b in unram) for a in inert)
        True

    There is no difference in providing a pAdicBase object and
    supplying a ring and prime ideal::

        sage: from modular_method.padics.pAdic_base import pAdicBase
        sage: from modular_method.padics.pAdic_extensions import ramified_quadratic_extensions
        sage: pAdics = pAdicBase(QQ, 2)
        sage: ramified_quadratic_extensions(ring=QQ, prime=2)
        [1, 3, 2, 6]
        sage: ramified_quadratic_extensions(pAdics=pAdics)
        [1, 3, 2, 6]

    """
    if pAdics is None:
        if prime is None:
            raise ValueError("At least the argument prime must be set")
        prime = Ideal(prime)
        if ring is None:
            ring = prime.ring()
        pAdics = pAdicBase(ring, prime)

    K = pAdics.number_field()
    result = [0, 0, 0, 0]
    result[0] = [K(1), pAdics.uniformizer()]
    # Representatives of <pi> / <pi^2>
    # which is all there is if char = 2 and
    # we are not interested in inert extensions
    if unramified and pAdics.characteristic() != 2:
        result[2] = [K(1)]
    if pAdics.characteristic() == 2 or inert or trivial:
        P = pAdics.prime_ideal()
        R.<x, y> = K[]
        T = pAdicTree([x, y], pAdics=pAdics)
        _, T._root = find_pAdic_roots(x*y, pAdics=pAdics,
                                      value_tree=T._root, precision=1)
        k = pAdics.valuation(2)
        if k > 0:
            Tunr, Tram = find_pAdic_roots(y - x^2, pAdics=pAdics,
                                          value_tree=T.root(), precision=2*k)
            Tunr = pAdicTree([x,y], root=Tunr)
            Tall = T.remove_variable('x')
            lsall = [N.representative()[0] for N in
                     Tall._root.children_at_level(2*k)]
            Tunr_ = Tunr.remove_variable('x')
            lsunr = [N.representative()[0] for N in
                     Tunr_._root.children_at_level(2*k)]
            if unramified:
                result[2] = lsunr
            lsquot = _quot_comp(lsall, lsunr, P^(2*k))
            lsrepr = [lsclass[0] for lsclass in lsquot]
            result[0] = [K(a1*a2) for a1 in result[0] for a2 in lsrepr]
        else:
            Tunr = T
        if inert or trivial:
            Ttriv, Tiner = find_pAdic_roots(y - x^2,
                                            pAdics=pAdics,
                                            value_tree=Tunr.root(),
                                            precision=2*k+1)
            Ttriv = pAdicTree([x, y], root=Ttriv)
            Ttriv_ = Ttriv.remove_variable('x')
            lstriv = [N.representative()[0] for N in
                      Ttriv_._root.children_at_level(2*k+1)]
        if trivial:
            result[3] = lstriv
        if inert:
            Tunr_ = Tunr.remove_variable('x')
            lsunr = [N.representative()[0] for N in
                     Tunr_._root.children_at_level(2*k+1)]
            lsquot = _quot_comp(lsunr, lstriv, P^(2*k+1))
            result[1] = [lsclass[0] for lsclass in lsquot]
    prresult = (True, inert, unramified, trivial)
    result = tuple(result[i] for i in range(4) if prresult[i])
    if len(result) == 1:
        return result[0]
    return result
