=====================================================================
 Computing the conductor exponent at 2 of :math:`E_{1, z, w}^\gamma`
=====================================================================

In this file we check that the conductors computed for :math:`E_{1, z,
w}^\gamma` as presented in Table~5.1 in the thesis are correct.

.. linkall

First we import the necessary methods from the `modular_method`
package.

::

   sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
   sage: from modular_method.elliptic_curves.tates_algorithm import tates_algorithm_multiple
   sage: from modular_method.padics.pAdic_base import pAdicBase

Next we compute the conductor exponent at `2` of :math:`E_{1, z,
w}^\gamma` for each choice of :math:`\gamma \in \\{1, -1, 2, -2\\}`
and store the results in the dictionary `N`. In the code below `E1` is
the curve :math:`E_{1, z, w}^\gamma` and `E2` is another curve to
which it is 2-isogenous. Using the method `tates_algorithm_multiple`
we can compute the conductor exponent of both these curves
simultaneously, and -- as the results are the same -- use the result
of the curve for which Tate's algorithm finishes first.

::

   sage: R.<z, w> = QQ[]
   sage: C = CoprimeCondition([z, w])
   sage: pAdics = pAdicBase(QQ, 2)
   sage: N = {}
   sage: for gamma in [1, -1, 2, -2]:
   ....:     E1 = EllipticCurve([0, 4*z*gamma, 0, 2*(z^2 + w)*gamma^2, 0])
   ....:     E2 = EllipticCurve([0, -8*z*gamma, 0, 8*(z^2 - w)*gamma^2, 0])
   ....:     N[gamma] = tates_algorithm_multiple(
   ....:         (E1, E2),
   ....:         pAdics=(pAdics, pAdics),
   ....:         initial_values=C.pAdic_tree(pAdics=pAdics),
   ....:         only_calculate=['conductor'],
   ....:         verbose=True,
   ....:         precision_cap=25,
   ....:     )

Note that the conditions stored in the dictionary `N` are not as
nicely formatted as the results in Table 5.1. To check they are the
same we need some additional methods. First of all a method that
iterates over all possible :math:`(z, w)` modulo an integer :math:`N`
that satisfy a condition `con`.

::

   sage: def elements_modulo(N, con):
   ....:     """Iterator over the z, w modulo `N` that satisfy `con`"""
   ....:     ls, N_ = con.pAdic_tree().give_as_congruence_condition()
   ....:     N_ = N_.gens_reduced()[0]
   ....:     dif = ceil(N / N_)
   ....:     for z, w in ls:
   ....:         for zdif in range(dif):
   ....:             for wdif in range(dif):
   ....:                 yield z + zdif*N_, w + wdif*N_

Next a method that given a list of tuples `(f, M)` with `f` a
polynomial and `M` a positive integer, and a condition `con`, tells us
if any of the :math:`(z, w)` satisfying `con` also satisfy `f(z, w) ==
0` modulo `M` for each tuple `(f, M)`.
   
::

   sage: def satisfy_congruences(congruences, con):
   ....:     """Give `True` iff any of the z, w satisfying `con` also satisfy the
   ....:     `congruences`

   ....:     The argument `congruences` is a list of tuples. Each tuple
   ....:     contains a polynomial $f$ and an integer $N$, and represents the
   ....:     congruence $f(z, w) \equiv 0$ (mod $N$).

   ....:     """
   ....:     maxN = lcm(N for poly, N in congruences)
   ....:     return any(
   ....:         all(mod(poly(z, w), N) == 0 for poly, N in congruences)
   ....:         for z, w in elements_modulo(maxN, con)
   ....:     )

A method which determines the possible values in a conditional value
for which at least one pair of :math:`(z, w)` that satisfies the
corresponding condition also satisfies `f(z, w) == 0` modulo `M` for
all tuples `(f, M)` in `congruences` with `f` a polynomial and `M` a
positive integer.
   
::

   sage: def possible_values(conditional_value, congruences):
   ....:     """Give the values of the `conditional_value` for which the
   ....:     `congruences` are satisfied

   ....:     """
   ....:     return [val for val, con in conditional_value if satisfy_congruences(congruences, con)]

Using the methods above we now construct a function that extracts for
each `\gamma` the possible values in `N` for which the corresponding
condition has at least one pair of :math:`(z, w)` satisfying both the
condition and `f(z, w) == 0` modulo `M` for all tuples `(f, M)` with
`f` a polynomial and `M` a positive integer in a given list
`congruences`.
   
::

   sage: def possible_N(congruences):
   ....:     """Give a dictionary of possible values satisfying the `congruences`
   ....:     for each twist

   ....:     """
   ....:     return {key: possible_values(value, congruences) for key, value in N.items()}

We now translate all conditions in Table 5.1 into congruences and
apply the method `possible_N` to show that for each of them and each
value of `\gamma` there is only one possible conductor
exponent. Furthermore those exponents agree with the table.

::

   sage: ({1: [8], -1: [8], 2: [8], -2: [8]} == possible_N([(w^2 - z^4 - 1, 2)]) and
   ....:  {1: [7], -1: [7], 2: [7], -2: [7]} == possible_N([(w^2 - z^4 - 8, 16)]) and
   ....:  {1: [4], -1: [3], 2: [6], -2: [6]} == possible_N([(w + z^2 - 8, 32), (z - 1, 4)]) and
   ....:  {1: [3], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 8, 32), (z - 3, 4)]) and
   ....:  {1: [2], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 24, 32), (z - 1, 4)]) and
   ....:  {1: [4], -1: [2], 2: [6], -2: [6]} == possible_N([(w + z^2 - 24, 32), (z - 3, 4)]) and
   ....:  {1: [6], -1: [6], 2: [4], -2: [2]} == possible_N([(w - z^2 - 8, 32), (z - 1, 4)]) and
   ....:  {1: [6], -1: [6], 2: [2], -2: [4]} == possible_N([(w - z^2 - 8, 32), (z - 3, 4)]) and
   ....:  {1: [6], -1: [6], 2: [3], -2: [4]} == possible_N([(w - z^2 - 24, 32), (z - 1, 4)]) and
   ....:  {1: [6], -1: [6], 2: [4], -2: [3]} == possible_N([(w - z^2 - 24, 32), (z - 3, 4)]) and
   ....:  {1: [5], -1: [5], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^4, 2^5)]) and
   ....:  {1: [6], -1: [6], 2: [5], -2: [5]} == possible_N([(w - z^2 - 2^4, 2^5)]) and
   ....:  {1: [3], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^5, 2^6), (z - 1, 4)]) and
   ....:  {1: [4], -1: [3], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^5, 2^6), (z - 3, 4)]) and
   ....:  {1: [6], -1: [6], 2: [4], -2: [3]} == possible_N([(w - z^2 - 2^5, 2^6), (z - 1, 4)]) and
   ....:  {1: [6], -1: [6], 2: [3], -2: [4]} == possible_N([(w - z^2 - 2^5, 2^6), (z - 3, 4)]) and
   ....:  {1: [3], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^6, 2^7), (z - 1, 4)]) and
   ....:  {1: [4], -1: [3], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^6, 2^7), (z - 3, 4)]) and
   ....:  {1: [6], -1: [6], 2: [4], -2: [3]} == possible_N([(w - z^2 - 2^6, 2^7), (z - 1, 4)]) and
   ....:  {1: [6], -1: [6], 2: [3], -2: [4]} == possible_N([(w - z^2 - 2^6, 2^7), (z - 3, 4)]) and
   ....:  {1: [0], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^7, 2^8), (z - 1, 4)]) and
   ....:  {1: [4], -1: [0], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^7, 2^8), (z - 3, 4)]) and
   ....:  {1: [6], -1: [6], 2: [4], -2: [0]} == possible_N([(w - z^2 - 2^7, 2^8), (z - 1, 4)]) and
   ....:  {1: [6], -1: [6], 2: [0], -2: [4]} == possible_N([(w - z^2 - 2^7, 2^8), (z - 3, 4)]) and
   ....:  {1: [1], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2, 2^8), (z - 1, 4)]) and
   ....:  {1: [4], -1: [1], 2: [6], -2: [6]} == possible_N([(w + z^2, 2^8), (z - 3, 4)]) and
   ....:  {1: [6], -1: [6], 2: [4], -2: [1]} == possible_N([(w - z^2, 2^8), (z - 1, 4)]) and
   ....:  {1: [6], -1: [6], 2: [1], -2: [4]} == possible_N([(w - z^2, 2^8), (z - 3, 4)]))
   True
