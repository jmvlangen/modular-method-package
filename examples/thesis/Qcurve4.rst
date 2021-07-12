=========================
 Code for Example 2.7.11
=========================

This file contains code that verifies the results from Example 2.7.11

.. linkall

The following import is required for the example to work

::

   sage: from modular_method import *

We start by constructing the Q-curve considered in the example.

::

   sage: K2.<sqrt2> = QuadraticField(2)
   sage: K5.<sqrt5> = QuadraticField(5)
   sage: from modular_method.number_fields.field_constructors import composite_field
   sage: K, phi2, phi5 = composite_field(K2, K5, give_maps=True)
   sage: sqrt2 = phi2(sqrt2)
   sage: sqrt5 = phi5(sqrt5)
   sage: sqrt10 = sqrt2 * sqrt5
   sage: a4 = -60*(15 + 10*sqrt2 + 5*sqrt5 + 2*sqrt10)
   sage: a6 = 80*(210 + 135*sqrt2 + 70*sqrt5 + 49*sqrt10)
   sage: E = Qcurve([0, 0, 0, a4, a6], guessed_degrees=[2, 3])

Next we name some elements of the Galois group: `s1` the trivial
element of the Galois group of `K` and `sn` the generator of the
Galois group of `K` over :math:`\QQ(\sqrt{n})`.

::

   sage: G = K.galois_group()
   sage: s1 = G(1)
   sage: s2 = next(s for s in G if s != s1 and s(sqrt2) == sqrt2)
   sage: s5 = next(s for s in G if s != s1 and s(sqrt5) == sqrt5)
   sage: s10 = next(s for s in G if s != s1 and s(sqrt10) == sqrt10)
   sage: Gls = [s1, s2, s5, s10]

Now we verify the degree map given in the example.

::

   sage: [E.degree_map(s) for s in Gls]
   [1, 2, 3, 6]

Next we verify that the splitting character is one of the characters
of conductor 15 and order 4.

::

   sage: D15 = [eps for eps in DirichletGroup(15) if eps.conductor() == 15 and eps.order() == 4]
   sage: len(D15)
   2
   sage: E.splitting_character() in D15
   True

We confirm that the splitting character field is indeed as
mentioned in the example.

::

   sage: L.<zeta15> = CyclotomicField(15)
   sage: Keps = L.subfield(zeta15 + zeta15^(-1))[0]
   sage: E.splitting_character_field().is_isomorphic(Keps)
   True
   sage: Keps.is_isomorphic(QQ[sqrt((15 - 3*sqrt(5))/2)])
   True

We check the splitting field is the same as the one mentioned in the
example.

::

   sage: Kbeta = composite_field(K, Keps)
   sage: E.splitting_field().is_isomorphic(Kbeta)
   True
   sage: Kbeta.is_isomorphic(QQ[sqrt(2), sqrt(15 - 3*sqrt(5))])
   True

Finally we check that the splitting image field is indeed as in the
example.

::

   sage: Lbeta = QQ[sqrt(sqrt(-1))*sqrt(2), sqrt(3)]
   sage: E.splitting_image_field().is_isomorphic(Lbeta)
   True
   sage: Lbeta.is_isomorphic(QQ[sqrt(-1), sqrt(3)])
   True
