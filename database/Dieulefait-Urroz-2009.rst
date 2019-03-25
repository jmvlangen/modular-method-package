============================================================================
Solving Fermat-type equations via modular Q-curves over polyquadratic fields
============================================================================

We run the computations for some results in the article "Solving
Fermat-type equations via modular Q-curves over polyquadratic fields"
written by Luis Dieulefait and Jorge Jiménez Urroz and published in
Journal für die reine und angewandte Mathematik, volume 633 (2009),
pages 183–195.

For this file the `arxiv version`_ of the article was used

.. _arxiv version: https://arxiv.org/abs/math/0611663
.. linkall

The article considers the equation :math:`x^4 + d y^2 = z^p` for which
we seek coprime integer solutions :math:`x`, :math:`y` and :math:`z`
with `p` a prime number and :math:`d` either 2 or 3.

Case d = 2
==========

In the case of the equation :math:`x^4 + 2 y^2 = z^p` the article
attaches to a primitive solution :math:`(A, B, C)` a Frey curve
defined over :math:`\QQ(\sqrt{-2})` which is also a Q-curve. We shall
denote ``Cp`` for the quantity :math:`C^p`.

::

   sage: load('load.sage')
   sage: R.<A, B> = QQ[]
   sage: Cp = A^4 + 2*B^2
   sage: con = CoprimeCondition([A, B]) & ~CongruenceCondition(A, 2)
   sage: K.<r> = QQ[sqrt(-2)]
   sage: a_invariants = [0, 4*A, 0, 2*(A^2 + r*B), 0]
   sage: G.<sigma> = K.galois_group()
   sage: isogenies = {sigma^0: (QQ(1), 1), sigma^1: (r, 2)}
   sage: E = FreyQcurve(a_invariants, isogenies=isogenies, condition=con)

We check the claims in the article that the field of complete
definition is unramified at 3 and that the curve has good or
semistable reduction at primes dividing 3.

::

   sage: 3.divides(E.complete_definition_field().discriminant())
   False
   sage: [E.has_additive_reduction(P) for P in K.primes_above(3)]
   [False, False]

We see as claimed in the article that all primes above odd prime
numbers are primes of semistable or good reduction.

::

   sage: Pbad = E.primes_of_possible_additive_reduction()
   Warning: Assuming that 160*A^2 + (-96*a)*B and 512*A^6 + (512*a)*A^4*B + 1024*A^2*B^2 + (1024*a)*B^3 are coprime outside ('(a)',).
   sage: [P.smallest_integer() for P in Pbad]
   [2]

As in the article we now distinguish two cases.

:math:`C` is not divisible by 3
-------------------------------

In this case we see that the restriction of scalars over the
decomposition field is an abelian surface of GL_2-type as mentioned in
the article.

::

   sage: E.does_decompose()
   True
   sage: E.splitting_image_field('conjugacy')
   (Number Field in a with defining polynomial x^2 - 2,)

We verify that the conductor exponent at the prime dividing 2 is
indeed 12 or 10. Furthermore the conductor of the restriction of
scalars and the resulting newform levels correspond with what we find
in the article.

::

   sage: P2 = K.prime_above(2)
   sage: E.conductor_exponent(P2)
   12 if ('A', 'B') == (1, 1) mod 2
   10 if ('A', 'B') == (1, 0) mod 2
   sage: E.conductor_restriction_of_scalars()
   2^(n0+6)*Norm(Rad_P( (512) * (A^2 + (-a)*B) * (A^2 + (a)*B)^2 ))
    where 
   n0 = 12 if ('A', 'B') == (1, 1) mod 2
        10 if ('A', 'B') == (1, 0) mod 2
   sage: E.newform_levels()
   [(512,)] if ('A', 'B') == (1, 1) mod 2
   [(256,)] if ('A', 'B') == (1, 0) mod 2

We compute the newforms at these levels and verify that there are
indeed 22 of them as the article states. Here we count all galois
conjugates hence for each conjugacy class we take the degree of the
corresponding coefficient field as the number of newforms in that
class.

::

   sage: nfs = E.newform_candidates()
   sage: sum(nf.coefficient_field().degree() for nf in nfs[0][0] + nfs[1][0])
   22

The article first eliminates all newforms with complex multiplication
under the assumption that :math:`p > 349`. As mentioned later on there
are only five remaining newforms up to conjugation all of level 512.

::

   sage: nfs = eliminate_cm_forms(E, nfs)
   sage: len(nfs[0][0]) + len(nfs[1][0])
   5
   sage: [nf[0].level() for nf in nfs[0][0] + nfs[1][0]]
   [512, 512, 512, 512, 512]

These newforms are now eliminated by noting that their third
coefficient should be congruent to some integer :math:`a_3` modulo a
prime above :math:`p` with :math:`|a_3| \le 2 \sqrt{3}`,
i.e. :math:`a_3 = 0, \pm 1, \pm 2, \pm 3`. By checking the primes in
the norm of the third coefficient of each newform minus each possible
:math:`a_3` we find that this can not be the case for :math:`p > 349`.

::

   sage: lcm(ZZ((nf[0].coefficient(3) - a3).absolute_norm()) for nf in nfs[0][0] + nfs[1][0]
   ....: for a3 in range(-3, 4)).prime_factors()
   [2, 3, 5, 7]

:math:`C` divisible by 3
------------------------

In this case only the way to determine whether the corresponding
galois representation is different. We thus get the same newforms.

::

   sage: nfs = E.newform_candidates()

Now for eliminating these forms the article suggests it is sufficient
to check that they all don't have a third coefficient that is
congruent to :math:`\pm 4` modulo a prime above :math:`p`. We check
this by computing all primes dividing the absolute norm of the third
coefficient of a newform minus :math:`\pm 4`.

::

   sage: lcm(ZZ((nf.coefficient(3) - a3).absolute_norm()) for nf in nfs[0][0] + nfs[1][0]
   ....: for a3 in [-4,4]).prime_factors()
   [2, 3, 5, 7, 17]

Case :math:`d = 3`
==================

The article next considers the equation :math:`x^4 + 3 y^2 = z^p`, for
which they use the same Frey curve only in this case with :math:`r` a
square root of -3. In this case for a primitive solution :math:`(A, B,
C)` we must assume that :math:`A` and :math:`B` are coprime and that
:math:`A` is not divisible by 3.

::

   sage: Cp = A^4 + 3*B^2
   sage: con = CoprimeCondition([A, B]) & ~CongruenceCondition(A, 3)
   sage: K.<r> = QuadraticField(-3)
   sage: a_invariants = [0, 4*A, 0, 2*(A^2 + r*B), 0]
   sage: G.<sigma> = K.galois_group()
   sage: L.<sqrtm2> = QuadraticField(-2)
   sage: isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm2, 2)}
   sage: E = FreyQcurve(a_invariants, isogenies=isogenies, condition=con)

The article reasons that the restriction of scalars of the curve
itself is not an abelian variety of GL_2-type which we verify.

::

   sage: E.does_decompose()
   False

The article then reasons that twisting the curve by :math:`\gamma =
2 + \sqrt{6}` would make it so the restriction of scalars of the curve
is an abelian variety of GL_2-type of dimension 4 which we verify. For
this we have to manually set the splitting character to the one
suggested in the article, the quadratic character of
:math:`\QQ(\sqrt{6})`.

::

   sage: gamma = 2 + QuadraticField(6).gen()
   sage: Ec = E.twist(gamma)
   sage: Ec._eps = {0 : [character_for_root(6)]}
   sage: Ec.does_decompose()
   True
   sage: [Lb.degree() for Lb in Ec.splitting_image_field('conjugacy')]
   [4]

After some computations and verifying modularity, the article states
that the possible levels of newforms corresponding to this curve 24,
96, 192 and 384 which we verify. We use the fact that the only bad
primes can be above 2 and 3.

::

   sage: Pbad = Ec.decomposition_field().primes_above(2*3)
   sage: Ec.newform_levels(bad_primes=Pbad)
   [(384,)] if ('A', 'B') == (0, 1) mod 2
   [(192,)] if ('A', 'B') == (1, 2), (3, 2) mod 4
   [(24,)]  if ('A', 'B') is 1 of 4 possibilities mod 8
   [(96,)]  if ('A', 'B') is 1 of 4 possibilities mod 8
   []       if ('A', 'B') == (1, 1) mod 2

As in the article we compute all newforms of these levels and first
eliminate all those newforms that have complex multiplication. We
check that the only newforms remaining are those of level 192 with
fifth coefficient squared equalt to 12 and those of level 384 with
seventh coefficient squared equal to -24 or -8.

::

   sage: nfs = Ec.newform_candidates(bad_primes=Pbad)
   sage: nfs = eliminate_cm_forms(Ec, nfs)
   sage: [nf[0].level() for nfsi in nfs for nf in nfsi[0]]
   [384, 384]
   sage: [nf[0].coefficient(7)^2 for nfsi in nfs for nf in nfsi[0]]
   [-8, -8]

We check as stated in the article that for each newform of level 384
the seventh coefficient is not congruent to :math:`z i` modulo primes
above :math:`p` for some integer :math:`z` of absolute value at
most 5. This we do by computing all prime numbers dividing the norm of
the different differences.

::

   sage: lcm(ZZ((nf[0].coefficient(3) - z*sqrt(nf[0].coefficient_field()(-1))).absolute_norm())
   ....: for nf in nfs[0][0] + nfs[1][0] for z in range(-5, 6)).prime_factors()
   [2, 3, 11, 19]

This shows the last result of the article.
