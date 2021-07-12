========================
 Code for Example 3.3.4
========================

This file contains the code for Example 3.3.4, so it can be verified
with SageMath's doctest system. There is also a `full example`_ for
the corresponding article.

.. _full example: ../literature/Bennett-Chen-2012.rst
.. linkall

The following import is required for the example to work

::

   sage: from modular_method import *

We first enter the conditions.

::

   sage: R.<a, b> = QQ[]
   sage: cl = a^2 + b^6
   sage: coprime = CoprimeCondition([a, b])
   sage: condition = coprime & PowerCondition(cl, 3)

Next we enter the Frey Q-curve.

::

   sage: K.<i> = QuadraticField(-1)
   sage: a_invariants = [0, 0, 0, -3*(5*b^3 + 4*a*i)*b,
   ....:                 2*(11*b^6 + 14*i*b^3*a - 2*a^2)]
   sage: E = FreyQcurve(a_invariants, condition=condition,
   ....:                guessed_degrees=[3])

We compute a decomposable twist.

::

   sage: Kdec = E.decomposition_field()
   sage: gamma = (-3 + sqrt(Kdec(-3))) / 2
   sage: E = E.twist(gamma)
   sage: E.does_decompose()
   True

Next we compute the corresponding newforms.

::

   sage: Kdef = E.definition_field()
   sage: nfs = E.newform_candidates(bad_primes=Kdef.primes_above(6))

We eliminate all those newforms with CM.

::

   sage: nfs = eliminate_cm_forms(E, nfs)

Introduce the second Frey curve and compute its corresponding
newforms.

::

   sage: E2 = FreyCurve([0, 0, 0, 3*b^2, 2*a], condition=condition)
   sage: nfs2 = E2.newform_candidates(bad_primes=[2, 3])

Eliminate newforms first for `E2` alone, then for both `E` and `E2`
together.

::

   sage: nfs2 = eliminate_by_traces(E2, nfs2, primes=[5, 7],
   ....:                            condition=coprime)
   sage: nfs_comb = combine_newforms(nfs, nfs2)
   sage: nfs_comb = eliminate_by_traces((E, E2), nfs_comb,
   ....:                                primes=[5, 7],
   ....:                                condition=coprime)

Eliminate the small primes and show nothing remains.

::

   sage: nfs_comb = eliminate_primes((E, E2), nfs_comb, 2*3*5*7)
   sage: nfs_comb
   []
