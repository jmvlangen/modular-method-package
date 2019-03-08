=========================================================================
 Multi-Frey Q-curve and the Diophantine equation :math:`a^2 + b^6 = c^n`
=========================================================================

We run the computations for some results in the article "Multi-Frey
Q-curve and the Diophantine equation :math:`a^2 + b^6 = c^n`" written
by Michael A. Bennett and Imin Chen and published in Algebra & Number
Theory, volume 6 (2012), no. 4. The article can be found at `Imin
Chen's website`_

.. _Imin Chen's website: http://people.math.sfu.ca/~ichen/pub/BeCh2.pdf
.. linkall

This article considers the equation :math:`a^2 + b^6 = c^l` for
coprime integers :math:`a` and :math:`b` and a prime number
:math:`l`. We will use the notation ``cl`` to denote the
:math:`l^{th}` power of :math:`c`. Furthermore the article assumes
that :math:`l` is at least 3.

::

   sage: R.<a,b> = QQ[]
   sage: cl = a^2 + b^6
   sage: C = (CoprimeCondition([a,b]) &
   ....:      PowerCondition(cl, 3))
   
The first Frey curve
====================

The authors associate to a solution of the equation a Frey curve
defined over the number field :math:`\QQ(i)` which is in fact a
Q-curve.

::

   sage: K.<i> = QuadraticField(-1)
   sage: L.<sqrtm3> = QuadraticField(-3)
   sage: G.<sigma> = K.galois_group()
   sage: a_invariants = [0, 0, 0, -3*(5*b^3 + 4*a*i)*b, 2*(11*b^6 + 14*i*b^3*a - 2*a^2)]
   sage: isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm3 , 3)}
   sage: E = FreyQcurve(a_invariants, isogenies=isogenies, condition=C)

According to the article the j-invariant of the curve is equal to
:math:`432 i \frac{b^3 (4 a - 5 i b^3)^3}{(a - i b^3) (a + i b^3)^3}`
      
::

   sage: E.j_invariant() == 432*i*(b^3*(4*a - 5*i*b^3)^3)/((a - i*b^3)*(a + i*b^3)^3)
   True

The discriminant should be :math:`-2^8 \cdot 3^3 (a - i b^3) (a + i
b^3)^3`.

::

   sage: E.discriminant() == -2^8*3^3*(a - i*b^3)*(a + i*b^3)^3
   True

The image of the degree map should be :math:`{1, 3}`.

::

   sage: E.degree_map_image() == [1, 3]
   True

The degree field should be :math:`\QQ(i)`.

::
   
   sage: E.degree_field().is_isomorphic(QQ[sqrt(-1)])
   True

A dual basis should be :math:`(-1, 3)`.

::

   sage: E.dual_basis() == ([-1], [3])
   True

The splitting character should be the product of the non trivial
character modulo 4 and the non trivial characte modulo 3.

::

   sage: (E.splitting_character() ==
   ....:  DirichletGroup(4).gen().extend(12) *
   ....:  DirichletGroup(3).gen().extend(12))
   True

The fixed field of the splitting character should be
:math:`\QQ(\sqrt{3})`

::

   sage: E.splitting_character_field().is_isomorphic(QQ[sqrt(3)])
   True

The fixed field of the splitting map should be :math:`\QQ(i,
\sqrt{3})`.

::

   sage: E.splitting_field().is_isomorphic(QQ[sqrt(3), sqrt(-1)])
   True

The image of the splitting map should lie in the same field.

::

   sage: E.splitting_image_field().is_isomorphic(QQ[sqrt(3), sqrt(-1)])
   True

In the article a twist of the curve is constructed to make sure the
restriction of scalars is the product of :math:`\QQ` simple varieties
of GL_2-type.

::

   sage: Kb = E.decomposition_field()
   sage: gamma = (-3 + sqrt(Kb(-3))) / 2
   sage: Eb = E.twist(gamma)
   sage: Eb.does_decompose()
   True
   sage: Eb.j_invariant() == E.decomposable_twist().j_invariant()
   True

It is shown in the article that the only primes of possible additive
reduction for the curve ``Eb`` are those above 2 and 3, for which
there is only one of each

::

   sage: q2 = Kb.prime_above(2)
   sage: q3 = Kb.prime_above(3)
   sage: Kb.primes_above(2*3) == [q2, q3]
   True

The conductor of ``Eb`` can be computed and agrees with the article.

::

   sage: N = Eb.conductor(additive_primes=[q2, q3]); N
   (4)*(1/4*izeta0^2 + 1)^n0*Rad_P( (-186624) * (b^3 + (-1/8*izeta0^3)*a) * (b^3 + (1/8*izeta0^3)*a)^3 )
    where 
   n0 = 0 if ('a', 'b') is 1 of 24 possibilities mod 9
        4 if ('a', 'b') is 1 of 48 possibilities mod 9

Also the conductor of its restriction of scalars agrees with the
article.

::

   sage: NR = Eb.conductor_restriction_of_scalars(additive_primes=[q2, q3]); NR
   65536*3^(2*n0+4)*Norm(Rad_P( (-186624) * (b^3 + (-1/28*izeta00zeta0^3 - 9/28*izeta00zeta0)*a) * (b^3 + (1/28*izeta00zeta0^3 + 9/28*izeta00zeta0)*a)^3 ))
    where 
   n0 = 0 if ('a', 'b') is 1 of 24 possibilities mod 9
        4 if ('a', 'b') is 1 of 48 possibilities mod 9

According to the article the restriction of scalars is itself a
:math:`\QQ` simple variety of GL_2-type.

::

   sage: Eb.number_of_splitting_maps(count_conjugates=False)
   1

Furthermore we can associate to it newforms of level 48 or 432.

::

   sage: Eb.newform_levels(bad_primes=[q2, q3])
   [(48,)]  if ('a', 'b') is 1 of 24 possibilities mod 9
   [(432,)] if ('a', 'b') is 1 of 48 possibilities mod 9

We get a list of newform candidates as presented in the article.

::

   sage: nfs = Eb.newform_candidates(bad_primes=[q2, q3], algorithm='magma')
   sage: F, = nfs[0][0]
   sage: G1, G2, G3 = nfs[1][0]
   sage: F.has_cm() and G1.has_cm() and G2.has_cm() and not G3.has_cm()
   True

The article has methods of eliminating the forms with complex
multiplication which we can do using a function

::

   sage: nfs = eliminate_cm_forms(Eb, nfs)
   sage: nfs[0][0] == []
   True
   sage: nfs[1][0] == [(G3, 0)]
   True

A second Frey curve
===================

In the article a second Frey curve is constructed to eliminate the
last newform ``G3``. This Frey curve is defined over the rationals.

::

   sage: a_invariants2 = [0, 0, 0, 3*b^2, 2*a]
   sage: E2 = FreyCurve(a_invariants2, condition=C)

As shown in the article this curve has discriminant :math:`-2^6 \cdot
3^3 (a^2 + b^6)`.

::

   sage: E2.discriminant() == (-2^6 * 3^3 * (a^2 + b^6))
   True

This curve can only have additive reduction at 2 and 3, hence we can
compute the conductor, which agrees with the result in the article.

::

   sage: N2 = E2.conductor(additive_primes=[2, 3]); N2
   2^n0*3^n1*Rad_P( (-1728) * (b^6 + a^2) )
    where 
   n0 =  6 if ('a', 'b') == (1, 0) mod 2
         5 if ('a', 'b') == (0, 1) mod 2
   n1 =  2 if ('a', 'b') is 1 of 24 possibilities mod 9
         3 if ('a', 'b') is 1 of 48 possibilities mod 9

Now we do some multi-Frey elimination on both curves using the primes
5 and 7 to compare traces at.

::

   sage: nfs2 = E2.newform_candidates(bad_primes=[2, 3], algorithm='magma')
   sage: nfs2 = eliminate_by_traces(E2, nfs2, primes=[5, 7], verbose=-1)
   sage: nfs_comb = combine_newforms(nfs, nfs2)
   sage: nfs_comb = eliminate_by_traces((Eb, E2), nfs_comb, primes=[5,7], verbose=-1)
   sage: lcm(f[2] for f in sum([ls[0] for ls in nfs_comb], [])).prime_factors()
   [2, 3, 5, 7]

This shows that the newform method eliminates all possible prime
exponents :math:`l` except for some small values which are discussed
seperately in the article.
