========================
 Code for Example 3.3.2
========================

This file contains the code for Example 3.3.2, so it can be verified
with SageMath's doctest system. There is also a `full example`_ for
the corresponding article.

.. _full example: ../literature/Bugeaud-Mignotte-Siksek-2008.rst
.. linkall

Entering the curves.

::

   sage: from modular_method import *
   sage: R.<psi> = QQ[]
   sage: con1 = (CongruenceCondition(psi - 8, 16) &
   sage:         CongruenceCondition(psi + 1, 5))
   sage: F1 = FreyCurve([0, 2*psi + 1, 0, psi^2 + psi, 0],
   ....:                condition=con1)
   sage: G2 = FreyCurve([0, 1, 0, -psi/4, 0],
   ....:                condition=con1)

Computing the conductor.

::

   sage: S = [2, 5]
   sage: F1.conductor(additive_primes=S)
   40*Rad_P( (16) * psi^2 * (psi + 1)^2 )
   sage: G2.conductor(additive_primes=S)
   160*Rad_P( (psi + 1) * psi^2 )

Computing the corresponding newforms.

::

   sage: nfs = [(f, g) for f in F1.newform_candidates(bad_primes=S)
   ....:        for g in G2.newform_candidates(bad_primes=S)]; nfs
   [(q + q^5 + O(q^6), q - 2*q^3 - q^5 + O(q^6)),
    (q + q^5 + O(q^6), q + 2*q^3 - q^5 + O(q^6)),
    (q + q^5 + O(q^6), q - a2*q^3 + q^5 + O(q^6))]

Eliminating newforms with Frobenius element at 3.

::

   sage: nfs = eliminate_by_trace((F1, G2), nfs, 3); nfs
   [(q + q^5 + O(q^6), q - 2*q^3 - q^5 + O(q^6), 0),
    (q + q^5 + O(q^6), q + 2*q^3 - q^5 + O(q^6), 12),
    (q + q^5 + O(q^6), q - a2*q^3 + q^5 + O(q^6), 12)]

Eliminating newforms with Frobenius element at 7.

::

   sage: nfs = eliminate_by_trace((F1, G2), nfs, 7); nfs
   [(q + q^5 + O(q^6), q - 2*q^3 - q^5 + O(q^6), 56),
    (q + q^5 + O(q^6), q + 2*q^3 - q^5 + O(q^6), 12),
    (q + q^5 + O(q^6), q - a2*q^3 + q^5 + O(q^6), 4)]

Second case of newforms.

::

   sage: nfs = [(f, g) for f in get_newforms(40)
   ....:        for g in get_newforms(20)]; nfs
   [(q + q^5 + O(q^6), q - 2*q^3 - q^5 + O(q^6))]

Elimination at multiple primes at once.

::

   sage: nfs = eliminate_by_traces((F1, G2), nfs, primes=[3, 7, 11, 13]); nfs
   [(q + q^5 + O(q^6), q - 2*q^3 - q^5 + O(q^6), 0)]
