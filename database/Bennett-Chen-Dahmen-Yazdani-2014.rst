=============================================
 On the equation :math:`a^3 + b^{3 n} = c^2`
=============================================

We run the computations for some results in the article "On the
equation :math:`a^3 + b^{3 n} = c^2`" written by M. A. Bennett, Imin
Chen, Sander R. Dahmen and Soroosh Yazdani and published in Acta
Arithmetica 163 (2014), no. 4, pages 327-343.

.. linkall

This article considers the equation :math:`a^2 + b^{3 l} = c^2` for
coprime integers :math:`a`, :math:`b` and :math:`c` and prime number
:math:`l`. We will use ``bl`` to denote the :math:`l^{th}` power of
:math:`b`.

::

   sage: load('load.sage')

First of all the article determines three different parametrizations
for the equation :math:`x^3 + y^3 = z^2` into coprime integers
:math:`s` and :math:`t`. With the help of these parametrizations the
article concludes that :math:`c` is not divisible by 3 and that there
remain three different valid parametrizations of :math:`b^l`.

::

   sage: R.<s, t> = QQ[]
   sage: coprime = CoprimeCondition([s, t])
   sage: bl1 = s^4 - 4*t*s^3 - 6*t^2*s^2 - 4*t^3*s + t^4
   sage: a1 = 2*(s^4 + 2*t*s^3 + 2*t^3*s + t^4)
   sage: c1 = 3*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)
   sage: C1 = coprime & ~CongruenceCondition(s - t, 2) & ~CongruenceCondition(s - t, 3)
   sage: bl2 = -3*s^4 + 6*t^2*s^2 + t^4
   sage: a2 = 3*s^4 + 6*t^2*s^2 - t^4
   sage: c2 = 6*s*t*(3*s^4 + t^4)
   sage: C2 = coprime & ~CongruenceCondition(s - t, 2) & ~CongruenceCondition(t, 3)
   sage: bl3 = 3*s^4 + 6*t^2*s^2 - t^4
   sage: a3 = -3*s^4 + 6*t^2*s^2 + t^4
   sage: c3 = 6*s*t*(3*s^4 + t^4)
   sage: C3 = coprime & ~CongruenceCondition(s - t, 2) & ~CongruenceCondition(t, 3)

Case :math:`c` odd
==================

The article starts with the case that :math:`c` is odd in which case
we have two Frey curves. The first being

::

   sage: E1 = FreyCurve([0, 0, 0, -3*a1, -2*c1], condition=C1)

We compute its invariants as in the article, which all agree with the
results in the article.

::

   sage: E1.discriminant() == -1728 * bl1^3
   True
   sage: E1.conductor()
   Warning: Assuming that s and t are coprime.
   64*3^n0*Rad_P( (-1728) * (s^4 - 4*s^3*t - 6*s^2*t^2 - 4*s*t^3 + t^4)^3 )
    where 
   n0 = 2 if ('s', 't') == (1, 2), (2, 1) mod 3
        3 if ('s', 't') is 1 of 4 possibilities mod 3

The second Frey curve is a Q-curve and is given by

::

   sage: K.<sqrt3> = QuadraticField(3)
   sage: a_invariants2 = [0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0]
   sage: E2 = FreyQcurve(a_invariants2, condition=C1, guessed_degrees=[2])

Again we compute the invariants as in the article, which agree with
the results in the article.

::

   sage: E2.discriminant() == (1664 - 960*sqrt3)*((s - t)^2 + 2*sqrt3*s*t)*((s - t)^2 - 2*sqrt3*s*t)^2
   True
   sage: G = K.galois_group()
   sage: [E2.splitting_map()(tau).minpoly() for tau in G]
   [x - 1, x^2 + 2]
   sage: E2.splitting_image_field().is_isomorphic(QQ[sqrt(-2)])
   True
   sage: (E2.splitting_character() ==
   ....:  DirichletGroup(4).gen().extend(12) *
   ....:  DirichletGroup(3).gen().extend(12))
   True
   sage: E2.conductor()
   Warning: Assuming that s and t are coprime.
   (64)*Rad_P( ((-960*sqrt3 + 1664)) * (s^2 + (2*sqrt3 - 2)*s*t + t^2) * (s^2 + (-2*sqrt3 - 2)*s*t + t^2)^2 )
   
The article computes the conductor of the galois representation
attached to the splitting map of this curve which is the square root
of the conductor of the restriction of scalars over the decomposition
field as this is in itself an abelian variety of GL_2 type.

::

   sage: E2.conductor_restriction_of_scalars()
   589824*Norm(Rad_P( ((-960*sqrt3 + 1664)) * (s^2 + (2*sqrt3 - 2)*s*t + t^2) * (s^2 + (-2*sqrt3 - 2)*s*t + t^2)^2 ))

We compute the newforms after level lowering, for which there are 10
conjugacy classes associated to the second curve according to the
article.

::

   sage: nfs1 = E1.newform_candidates(bad_primes=[2, 3], algorithm='magma')
   sage: nfs2 = E2.newform_candidates(bad_primes=K.primes_above(2*3), algorithm='magma')
   sage: len(nfs2)
   10

We apply the multi-Frey method as in the article. First we do a trick
to put the newform lists in the right format.

::

   sage: nfs1e = eliminate_by_traces(E1, nfs1, primes=[])
   sage: nfs2e = eliminate_by_traces(E2, nfs2, primes=[])
   sage: nfs = combine_newforms(nfs1e, nfs2e)
   sage: nfs = eliminate_by_traces((E1, E2), nfs, primes=[5, 7, 11])

According to the article all the remaining cases for primes :math:`l`
smaller than 13 are in an explicit list. We check that this is true.

::

   sage: nfs = eliminate_primes((E1, E2), nfs, 2*3*5*7*11)
   sage: F1, F2, F3, F4, F5, F6, F7, F8, F9, F10 = nfs2
   sage: g11 = nfs1[0][0][5]
   sage: g12 = nfs1[1][0][17]
   sage: g13 = nfs1[1][0][26]
   sage: nfs[0][0] == [(g11, F1, 0), (g11, F2, 0), (g11, F4, 0), (g11, F5, 0)]
   True
   sage: nfs[1][0] == [(g12, F3, 0), (g12, F6, 0), (g13, F3, 0), (g13, F6, 0)]
   True

These newforms are eliminated in the article using an image of inertia
argument and the image of the projectivized galois representation
respectively.

Case :math:`c` even
===================

In this case the article uses three Frey-Hellegouarch curves. Note
that all of these in fact have two different choices, corresponding to
the choice of parametrization of :math:`b^l`. We start with the first
curve.

::

   sage: E11 = FreyCurve([0, 0, 0, -3*a2, -2*c2], condition=C2)
   sage: E12 = FreyCurve([0, 0, 0, -12*a3, -16*c3], condition=C3)

These curves have the same conductor according to the article and are
precisely as we compute here.

::

   sage: E11.conductor()
   Warning: Assuming that s and t are coprime.
   32*3^n0*Rad_P( (1728) * (3*s^4 - 6*s^2*t^2 - t^4)^3 )
    where 
   n0 = 2 if ('s', 't') == (0, 1), (0, 2) mod 3
        3 if ('s', 't') is 1 of 4 possibilities mod 3
   sage: E12.conductor()
   Warning: Assuming that s and t are coprime.
   32*3^n0*Rad_P( (-110592) * (3*s^4 + 6*s^2*t^2 - t^4)^3 )
    where 
   n0 = 2 if ('s', 't') == (0, 1), (0, 2) mod 3
        3 if ('s', 't') is 1 of 4 possibilities mod 3

The second Frey curve introduced in this case is

::

   sage: K.<sqrt3> = QuadraticField(3)
   sage: a_invariants21 = [0, 4*(sqrt3 - 1)*t, 0, -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 - sqrt3)*t^2), 0]
   sage: a_invariants22 = [0, 4*(sqrt3 - 1)*t, 0, -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 + sqrt3)*t^2), 0]
   sage: E21 = FreyQcurve(a_invariants21, condition=C2, guessed_degrees=[2])
   sage: E22 = FreyQcurve(a_invariants22, condition=C3, guessed_degrees=[2])

Again both curves have the same conductor and the conductor according
to the article is the same as computed here.

::

   sage: E21.conductor()
   Warning: Assuming that s and t are coprime.
   (64)*Rad_P( ((39936*sqrt3 - 69120)) * (s^2 + (2/3*sqrt3 - 1)*t^2) * (s^2 + (-2/3*sqrt3 - 1)*t^2)^2 )
   sage: E22.conductor()
   Warning: Assuming that s and t are coprime.
   (64)*Rad_P( ((39936*sqrt3 - 69120)) * (s^2 + (2/3*sqrt3 + 1)*t^2) * (s^2 + (-2/3*sqrt3 + 1)*t^2)^2 )

Furthermore both curve have a restriction of scalar that is an abelian
variety of GL_2-type. In the article they compute the conductor of a
galois representation attached to a splitting map, which is again the
square root of the conductor of this restriction of scalar. This
agrees with the following computation.

::

   sage: E21.conductor_restriction_of_scalars()
   589824*Norm(Rad_P( ((39936*sqrt3 - 69120)) * (s^2 + (2/3*sqrt3 - 1)*t^2) * (s^2 + (-2/3*sqrt3 - 1)*t^2)^2 ))
   sage: E22.conductor_restriction_of_scalars()
   589824*Norm(Rad_P( ((39936*sqrt3 - 69120)) * (s^2 + (2/3*sqrt3 + 1)*t^2) * (s^2 + (-2/3*sqrt3 + 1)*t^2)^2 ))

As in the article we now apply the multi-Frey method to these two/four
Frey curves first.

::

   sage: nfs11 = E11.newform_candidates(bad_primes=[2,3], algorithm='magma')
   sage: nfs12 = E12.newform_candidates(bad_primes=[2,3], algorithm='magma')
   sage: nfs21 = E21.newform_candidates(bad_primes=K.primes_above(2*3), algorithm='magma')
   sage: nfs22 = E22.newform_candidates(bad_primes=K.primes_above(2*3), algorithm='magma')
   sage: nfs11e = eliminate_by_traces(E11, nfs11, primes=[])
   sage: nfs12e = eliminate_by_traces(E12, nfs12, primes=[])
   sage: nfs21e = eliminate_by_traces(E21, nfs21, primes=[])
   sage: nfs22e = eliminate_by_traces(E22, nfs22, primes=[])
   sage: nfsc1 = combine_newforms(nfs11e, nfs21e)
   sage: nfsc2 = combine_newforms(nfs12e, nfs22e)
   sage: nfsc1 = eliminate_by_traces((E11, E21), nfsc1, primes=[5, 7, 11])
   sage: nfsc2 = eliminate_by_traces((E12, E22), nfsc2, primes=[5, 7, 11])

In the article they only consider remaining cases for :math:`l \ge
13`, hence we eliminate all other cases and find the same remaining
cases as in the article.

::

   sage: nfsc1 = eliminate_primes((E11, E21), nfsc1, 2*3*5*7*11)
   sage: nfsc2 = eliminate_primes((E12, E22), nfsc2, 2*3*5*7*11)
   sage: F1, F2, F3, F4, F5, F6, F7, F8, F9, F10 = nfs21
   sage: g11 = nfs11[0][0][0]
   sage: g12 = nfs11[1][0][1]
   sage: g13 = nfs11[1][0][5]
   sage: nfsc1[0][0] == [(g11, F1, 0), (g11, F2, 0), (g11, F4, 0), (g11, F5, 0)]
   True
   sage: nfsc1[1][0] == [(g12, F3, 0), (g12, F6, 0), (g13, F3, 0), (g13, F6, 0)]
   True
   sage: F1, F2, F3, F4, F5, F6, F7, F8, F9, F10 = nfs22
   sage: g11 = nfs12[0][0][0]
   sage: nfsc2[0][0] == [(g11, F1, 0), (g11, F2, 0), (g11, F4, 0), (g11, F5, 0)]
   True
   sage: nfsc2[1][0] == []
   True

The article introduces a third pair of curves to eliminate some
remaining newforms.

::

   sage: a_invariants31 = [0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3+3)*s^2), 0]
   sage: a_invariants32 = [0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3-3)*s^2), 0]
   sage: E31 = FreyQcurve(a_invariants31, condition=C2, guessed_degrees=[2])
   sage: E32 = FreyQcurve(a_invariants32, condition=C3, guessed_degrees=[2])

We compute the newforms of these curves and quickly note that their
level is indeed the indicated level in the article.

::

   sage: nfs31 = E31.newform_candidates(bad_primes=K.primes_above(2*3), algorithm='magma')
   sage: nfs32 = E32.newform_candidates(bad_primes=K.primes_above(2*3), algorithm='magma')
   sage: nfs31[0].level()
   2304
   sage: nfs32[0].level()
   2304

Next we perform the multi-Frey method as indicated in the article and
check we indeed get the same cases as indicated.

::

   sage: nfs31e = eliminate_by_traces(E31, nfs31, primes=[])
   sage: nfsc1 = combine_newforms(nfsc1, nfs31e)
   sage: nfsc1 = eliminate_by_traces((E11, E21, E31), nfsc1, primes=[5, 7, 11])
   sage: nfsc1 = eliminate_primes((E11, E21, E31), nfsc1, 2*3*5*7*11)
   sage: F1, F2, F3, F4, F5, F6, F7, F8, F9, F10 = nfs21
   sage: G1, G2, G3, G4, G5, G6, G7, G8, G9, G10 = nfs31
   sage: g12 = nfs11[1][0][1]
   sage: g13 = nfs11[1][0][5]
   sage: nfsc1[1][0] == [(g12, F3, G5, 0), (g12, F3, G6, 0), (g12, F6, G7, 0), (g12, F6, G8, 0), (g13, F3,
   ....:  G7, 0), (g13, F3, G8, 0), (g13, F6, G5, 0), (g13, F6, G6, 0)]
   True

The rest of the cases is now treated separately by the article.
