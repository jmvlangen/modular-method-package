==============================================================
 Winding quotients and some variants of Fermat's Last Theorem
==============================================================

We run some of the computations done in the article "Winding quotients
and some variants of Fermat's Last Theorem" written by Henri Darmon
and Loïc Merel and published in Journal für die reine und angewandte
Mathematik, volume 490 (1997), pages 81-100.

The article can be found on the `website of Merel`_

.. _website of Merel: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
.. linkall

This article considers integer solutions :math:`(x, y, z)` to the
equations

1) :math:`x^n + y^n = 2 z^n`,
2) :math:`x^n + y^n = z^2`, and
3) :math:`x^n + y^n = z^3`,

with a strictly positive integer :math:`n`. The article only considers
primitive solutions, i.e. :math:`gcd(x, y, z) = 1`, which are
non-trivial, i.e. :math:`|x y z| > 1`.

We note that the article in the introduction already gives examples of
why small cases for :math:`n` have infinitely many solutions. The
non-existence of non-trivial primitive solutions will be proven for
:math:`n \ge 3` and :math:`n \ge 4` in the case of equation 2. The
article however remarks that the method discussed in the article will
only have to work for :math:`n = p \ge 7` prime, due to previous
results.

The article denotes a non-trivial primite solution to any of the three
equations by :math:`(a, b, c)`. We will also introduce the notation
``ap``, ``bp`` and ``cp`` to denote the respective variable to the
:math:`p^{th}` power. The variable ``c0`` is an auxiliary variable the
article introduces for one of the Frey curves for equation 3.

::

   sage: load('load.sage')
   sage: R.<a, b, c, ap, bp, cp, c0> = QQ[]
   
For equation 1 the article introduces a single Frey curve. Furthermore
the article assumes without loss of generality that :math:`a` is -1
modulo 4.

::

   sage: S1 = QQ[ap, cp]
   sage: C1 = (CoprimeCondition([ap, cp]) & CongruenceCondition(ap + 1, 4) &
   ....:       PowerCondition(ap, 4) & PowerCondition(2*cp - ap, 4) & PowerCondition(cp, 4))
   sage: E1 = FreyCurve([0, S1(-ap - 2*cp), 0, S1((-ap)*(-2*cp)), 0], condition=C1)

For equation 2 the article introduces 2 different Frey curves, one for
the case that :math:`a b` is even and one for the case it is odd. In
the first case the article assume without loss of generality that
:maht:`a` is even and :math:`c \equiv 1` modulo 4. In the second case
the article assumes :math:`a \equiv -1` modulo 4. For the curve `E21`
we will remove a power of 2 from `ap` to reduce on computation time
for the conductor.

::

   sage: S2 = QQ[ap, c]
   sage: ap_ = 2^6 * ap
   sage: C21 = (CoprimeCondition([ap_, c]) & CongruenceCondition(c - 1, 4) &
   ....:        PowerCondition(ap_, 7) & PowerCondition(c^2 - ap_, 2))
   sage: E21 = FreyCurve([1, S2((c - 1)/4), 0, S2(ap), 0], condition=C21)
   sage: C22 = (CoprimeCondition([ap, c]) & ~CongruenceCondition(ap*(c^2 - ap), 2) &
   ....:        CongruenceCondition(ap + 1, 4) & PowerCondition(ap, 2) &
   ....:        PowerCondition(c^2 - ap, 2))
   sage: E22 = FreyCurve([0, S2(2*c), 0, S2(ap), 0], condition=C22)

For equation 3 the equation again introduces two Frey curves depending
on the parity of :math:`a b`. In case :math:`c` is even the article
writes :math:`c = 2 c_0`. Otherwise it assumes without loss of
generality that :math:`a` is odd and :math:`b` is even.

::

   sage: S3 = QQ[bp, c]
   sage: C31 = (CoprimeCondition([bp, c]) & ~CongruenceCondition(bp*(c^3 - bp), 2) &
   ....:        CongruenceCondition(c, 2) & PowerCondition(bp, 4) &
   ....:        PowerCondition((2*c)^3 - bp, 4))
   sage: E31 = FreyCurve([0, 0, S3(bp), S3(-3*((c/2)^3 + bp)*(c/2)),
   ....:                  S3(-(c/2)^3*(2*(c/2)^3 - 5*bp))],
   ....:                 condition=C31)
   sage: C32 = (CoprimeCondition([bp, c]) & CongruenceCondition(bp, 2) &
   ....:        PowerCondition(bp, 4) & PowerCondition(c^3 - bp, 4))
   sage: E32 = FreyCurve([S3(c), S3(-c^2), 0, S3(-3/2*c*bp), S3(bp*(c^3 + bp/4))],
   ....:                 condition=C32)

We check that the discriminants are indeed as listed in the article.

::

   sage: E1.discriminant() == 2^6*(ap*(2*cp - ap)*cp)^2
   True
   sage: E21.discriminant() == ap^2*(c^2 - 2^6*ap)
   True
   sage: E22.discriminant() == 2^6*ap^2*(c^2 - ap)
   True
   sage: E31.discriminant() == 3^3*((2*c0)^3 - bp)^3*bp
   True
   sage: E32.discriminant() == 3^3*(c^3 - bp)^3*bp
   True

In proposition 1.1 the conductors of the different curves are
given. We check that these are correct and see that in fact the
article misses the case for the last two curves where the conductor
exponent at 3 may be 2. This can be verified against the table of
Papadopoulus to exist in case 3 divides :math:`b` or the property
:math:`P_2` holds.

::

   sage: E1.conductor(additive_primes=[2])
   2^n0*Rad_P( (64) * cp^2 * ap^2 * (ap - 2*cp)^2 )
    where 
   n0 = 5 if ('ap', 'cp') == (3, 1), (3, 3) mod 4
        1 if ('ap', 'cp') is 1 of 32 possibilities mod 128
   sage: E21.conductor(additive_primes=[2])
   2*Rad_P( ap^2 * (c^2 - 64*ap) )
   sage: E22.conductor(additive_primes=[2])
   32*Rad_P( (64) * ap^2 * (c^2 - ap) )
   sage: E31.conductor(additive_primes=[3])
   3^n0*Rad_P( (27) * bp * (8*c0^3 - bp)^3 )
    where 
   n0 = 2 if ('bp', 'c0') is 1 of 710046 possibilities mod 2187
        3 if ('bp', 'c0') is 1 of 24 possibilities mod 9
        1 if ('bp', 'c0') is 1 of 1458 possibilities mod 2187
   sage: E32.conductor(additive_primes=[3])
   3^n0*Rad_P( (27) * bp * (c^3 - bp)^3 )
    where 
   n0 = 2 if ('bp', 'c') is 1 of 710046 possibilities mod 2187
        3 if ('bp', 'c') is 1 of 24 possibilities mod 9
        1 if ('bp', 'c') is 1 of 1458 possibilities mod 2187

We verify collary 3.2 by computing the newforms associated to the
different curves.

::

   sage: nfs1 = E1.newform_candidates(bad_primes=[2]); nfs1
   [q - 2*q^5 + O(q^6)] if ('ap', 'cp') == (3, 1), (3, 3) mod 4
   []                   if ('ap', 'cp') is 1 of 32 possibilities mod 128
   sage: nfs21 = E21.newform_candidates(bad_primes=[2]); nfs21
   []
   sage: nfs22 = E22.newform_candidates(bad_primes=[2]); nfs22
   [q - 2*q^5 + O(q^6)]
   sage: nfs31 = E31.newform_candidates(bad_primes=[3]); nfs31
   []                   if ('bp', 'c0') is 1 of 711504 possibilities mod 2187
   [q - 2*q^4 + O(q^6)] if ('bp', 'c0') is 1 of 24 possibilities mod 9
   sage: nfs32 = E32.newform_candidates(bad_primes=[3]); nfs32
   []                   if ('bp', 'c') is 1 of 711504 possibilities mod 2187
   [q - 2*q^4 + O(q^6)] if ('bp', 'c') is 1 of 24 possibilities mod 9

The last thing we check is that the remaining newforms have complex
multiplication as claimed in the article.

::

   sage: nfs1[0][0][0].has_cm()
   True
   sage: nfs22[0].has_cm()
   True
   sage: nfs31[1][0][0].has_cm()
   True
   sage: nfs32[1][0][0].has_cm()
   True
