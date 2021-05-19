========================
 Code for Example 2.5.7
========================

This file contains the code for the computations in Example 2.5.7.

.. linkall

We start by introducing the Q-curve.

::

   sage: _.<sqrt17> = QuadraticField(17)
   sage: E = Qcurve([0, 12, 0, 18*(1 + sqrt17), 0],
   ....:            guessed_degrees=[2])

Next we do the computations in the example to determine the complete
definition field and the matrix of :math:`c_E c_\beta^{-1}`.

::

   sage: K = E.complete_definition_field()
   sage: sqrtm2, sqrt17 = sqrt(K(-2)), sqrt(K(17))
   sage: K.is_isomorphic(QQ[sqrtm2, sqrt17])
   True
   sage: G = K.galois_group()
   sage: s2 = next(s for s in G if s != G(1) and s(sqrtm2) == sqrtm2)
   sage: s17 = next(s for s in G if s != G(1) and s(sqrt17) == sqrt17)
   sage: matrix([[E.c(s, t) / E.c_splitting_map(s, t)
   ....:          for t in [G(1), s2, s17, s2*s17]]
   ....:         for s in [G(1), s2, s17, s2*s17]])
   Warning: The restriction of scalars of this Q-curve over the decomposition field does not decompose into abelian varieties of GL_2-type. Use the method decomposable_twist to find a twist that does.
   [ 1  1  1  1]
   [ 1 -1  1 -1]
   [ 1 -1  1 -1]
   [ 1  1  1  1]

We compute the generators of the unit group.

::

   sage: u0, u1 = K.unit_group().gens_values()
   sage: u0 == -1 and u1^(-1) == 4 + sqrt17
   True
   sage: u0, u1 = u0, u1^(-1)

TODO: verify the inconsistent system of equations

We verify the computation of :math:`S` in the example.

::

   sage: S = E._decomposable_twist_set()
   sage: S.reverse()
   sage: S == K.primes_above(2)
   True

Finally we repeat the last computations to verify the remaining
results in the example.

::

   sage: gamma = E._decomposable_twist()
   sage: alpha = {s : sqrt(s(gamma) / gamma) for s in G}
   sage: (alpha[G(1)] == 1 and
   ....:  alpha[s2] == (7*sqrt17 - 29) / (2*sqrtm2) and
   ....:  alpha[s17] == 1 and
   ....:  alpha[s2*s17] == (7 * sqrt17 - 29) / (2*sqrtm2))
   True
   sage: all(P in S for P, _ in K.ideal(alpha[s2]).factor())
   True
   sage: all(E.c(s, t) / E.c_splitting_map(s, t) ==
   ....:     alpha[s] * s(alpha[t]) / alpha[s*t]
   ....:     for s in G for t in G)
   True
