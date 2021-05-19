========================
 Code for Example 3.3.3
========================

This file contains the code for Example 3.3.3, so it can be verified
with SageMath's doctest system. There is also a `full example`_ for
the corresponding article.

.. _full example: ../literature/Kraus-1998.rst
.. linkall

We enter the condition.

::

   sage: from modular_method import *
   sage: R.<a, b> = QQ[]
   sage: cl = a^3 + b^3
   sage: coprime = CoprimeCondition([a, b])
   sage: condition = (coprime &
   ....:              CongruenceCondition(a - 2, 4) &
   ....:              CongruenceCondition(b - 1, 4) &
   ....:              CongruenceCondition(cl, 3) &
   ....:              PowerCondition(cl, 5))

We construct the Frey curve and compute its newforms.

::

   sage: Eab = FreyCurve([0, 0, 0, 3*a*b, b^3 - a^3], condition=condition)
   sage: nfs = Eab.newform_candidates(); nfs
   Warning: Assuming that a and b are coprime.
   Warning: The bad primes chosen by default only take into account primes of additive reduction.
   [q + 2*q^5 + O(q^6)]

We construct the polynomials that are l-th powers.

::

   sage: K.<zeta3> = CyclotomicField(3)
   sage: poly0 = 3*(a + b)
   sage: poly1 = (a + zeta3*b) / (1 - zeta3)
   sage: poly2 = (a + zeta3^2*b) / (1 - zeta3)

We perform the Kraus method.

::

   sage: kraus_method(Eab, nfs, 5, (poly0, poly1, poly2),
   ....:              primes=prime_range(7, 50), condition=coprime)
   [(q + 2*q^5 + O(q^6), 0)]

Now again after a little modification so it works.

::

   sage: nfs5 = [(nf, 5) for nf in nfs]
   sage: kraus_method(Eab, nfs5, 5, (poly0, poly1, poly2),
   ....:              primes=prime_range(7, 50), condition=coprime)
   []
