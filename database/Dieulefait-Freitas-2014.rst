==============================================================================================
 The Fermat-type Equations :math:`x^5 + y^5 = 2 z^p` or :math:`3 z^p` solved through Q-curves
==============================================================================================

We run the computations for some results in the article "The
Fermat-type Equations :math:`x^5 + y^5 = 2 z^p` or :math:`3 z^p`
solved through Q-curves" written by Luis Dieulefait and Nuno Freitas
and published in Mathematics of Computation, volume 83 (2014), no.
286, pages 917-933.

.. linkall

This article considers integer solutions to equations of the form
:math:`x^5 + y^5 = d z^l` for :math:`d` either 2 or 3 and math:`l` a
prime number. Solutions :math:`a, b, c` to this equation will be
integers with :math:`a` and :math:`b` coprime. We will write ``dcl``
for the quantity :math:`d c^l`

::

   sage: load('load.sage')
   sage: R.<a, b> = QQ[]
   sage: dcl = a^5 + b^5
   sage: C = CoprimeCondition([a, b])

Furthermore, the article makes use of the factorization of
:math:`x^5 + y^5` over the field :math:`\QQ(\sqrt{5})`. In particular
writing ``w`` for the root of the polynomial :math:`x^2 + x - 1` it
uses the following factors.

::

   sage: S.<x> = QQ[]
   sage: K.<w> = NumberField(x^2 + x - 1)
   sage: G.<sigma> = K.galois_group()
   sage: phi1 = a^2 + w*a*b + b^2
   sage: phi2 = phi1.change_ring(sigma.as_hom())
   sage: dcl == (a+b) * phi1 * phi2
   True

From these factors the article constructs a Frey curve that is also a
Q-curve.

::

   sage: a_invariants = [0, 2*(a + b), 0, -sigma(w)*phi1, 0]
   sage: L.<sqrtm2> = QuadraticField(-2)
   sage: isogenies = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm2, 2)}
   sage: E = FreyQcurve(a_invariants, isogenies=isogenies, condition=C)

We check that the discriminant is indeed as stated in the article.

::

   sage: E.discriminant() == 2^6 * sigma(w) * phi1^2 * phi2
   True

Most computations on the Q-curve itself are done in section 6 of the
article. We check here that these results agree with our computation.

::

   sage: E.degree_field().is_isomorphic(QQ[sqrt(5)])
   True
   sage: E.dual_basis() == ([5], [2])
   True

As in the article we will call the field of complete definition of the
curve :math:`L` and define the elements :math:`\sigma` and
:math:`\tau` of its galois group in a similar way.

::

   sage: L = E.complete_definition_field()
   sage: L.is_isomorphic(QQ[sqrt(5), sqrt(-2)])
   True
   sage: G.<sigma, tau> = L.galois_group()
   sage: sigma(sqrt(L(5))) == sqrt(L(5))
   True
   sage: sigma(sqrt(L(-2))) == -sqrt(L(-2))
   True
   sage: tau(sqrt(L(5))) == -sqrt(L(5))
   True
   sage: tau(sqrt(L(-2))) == sqrt(L(-2))
   True

Next we compute the values of the associated cocycle, which agrees
with table 2 in the article.

::

   sage: S = [G(1), sigma, tau, sigma*tau]
   sage: matrix([[E.c(g, h) for h in S] for g in S])
   [ 1  1  1  1]
   [ 1  1 -1 -1]
   [ 1  1 -2 -2]
   [ 1  1  2  2]

Our computations give the same result for :math:`\xi(C)_{\pm}` and the
splitting character.

::

   sage: E.xi_pm()
   [(5, 2)]
   sage: (E.splitting_character() ==
   ....:  DirichletGroup(4).gen().extend(20) *
   ....:  DirichletGroup(5).gen().extend(20))
   True

The fixed field of the splitting character is :math:`\QQ` adjoint
:math:`\theta = \sqrt{\frac{1}{2}(5 + \sqrt{5})}`, which is a root of
:math:`x^4 - 5 x^2 + 5`. According to the article this is also the
fixed field of the splitting map.

::

   sage: Ke.<theta> = NumberField(x^4 - 5*x^2 + 5)
   sage: E.splitting_character_field().is_isomorphic(Ke)
   True
   sage: E.splitting_field().is_isomorphic(Ke)
   True

According to the article this curve does not decompose over the
decomposition field, but if we twist with the element :math:`\gamma =
2 \theta^2 - \theta - 5` it does. The latter even decomposes over the
field :math:`\QQ(\theta)`. Furthermore it decomposes as the the
product of two non-isogenous abelian surfaces of GL_2-type.

::

   sage: E.does_decompose()
   False
   sage: gamma = 2*theta^2 - theta - 5
   sage: Ec = E.twist(gamma)
   sage: Ec.does_decompose()
   True
   sage: Ec.decomposition_field().is_isomorphic(Ke)
   True
   sage: Ec.number_of_splitting_maps(count_conjugates=False)
   2
   sage: Ec.splitting_image_field('conjugacy')
   (Number Field in zeta80 with defining polynomial x^2 + 2*x + 2,
    Number Field in zeta80 with defining polynomial x^2 + 2*x + 2)

Now we again check that ``Ec`` has the invariants as mentioned in
section 3.1 of the article. Note that the invariant :math:`c_4` as
printed in the article is wrong, as the second - should be a +.

::

   sage: iota = K.embeddings(Ke)[0]
   sage: iso = Ec.definition_field().embeddings(Ke)[0]
   sage: bar = K.galois_group().gen()
   sage: Ec.discriminant().change_ring(iso) == gamma^6 * 2^6 * (bar(w) * phi1^2 * phi2).change_ring(iota)
   True
   sage: Ec.c4().change_ring(iso) == -gamma^2 * 2^4 * (bar(w)*phi1 + 2^2*w*phi2).change_ring(iota)
   True
   sage: Ec.c6().change_ring(iso) == -gamma^3 * 2^6 * (a+b)*(bar(w)*phi1 - 2^3*w*phi2).change_ring(iota)
   True

As in the article we denote the only primes above 2 and 5 by ``B2``
and ``B5`` respectively. The conductor exponent at ``B5`` we compute
is the same as presented in proposition 3.4.

::

   sage: B2 = Ec.definition_field().prime_above(2)
   sage: B5 = Ec.definition_field().prime_above(5)
   sage: Ec.conductor_exponent(B5)
   2 if ('a', 'b') is 1 of 20 possibilities mod 5
   0 if ('a', 'b') is 1 of 4 possibilities mod 5
   sage: Ec.conductor_exponent(B5)[1][1]
   The condition that ('a', 'b') == (1, 4), (2, 3), (3, 2), (4, 1) mod 5

The conductor exponent at ``B2`` is the same as presented in
proposition 3.5.

::

   sage: Ec.conductor_exponent(B2)
   8 if ('a', 'b') is 1 of 6 possibilities mod 4
   6 if ('a', 'b') is 1 of 4 possibilities mod 4
   4 if ('a', 'b') is 1 of 4 possibilities mod 8
   0 if ('a', 'b') is 1 of 4 possibilities mod 8
   sage: Ec.conductor_exponent(B2)[0][1]
   The condition that ('a', 'b') == (1, 1), (1, 2), (2, 1), (2, 3), (3, 2), (3, 3) mod 4
   sage: Ec.conductor_exponent(B2)[1][1]
   The condition that ('a', 'b') == (0, 1), (0, 3), (1, 0), (3, 0) mod 4
   sage: Ec.conductor_exponent(B2)[2][1]
   The condition that ('a', 'b') == (1, 7), (3, 5), (5, 3), (7, 1) mod 8
   sage: Ec.conductor_exponent(B2)[3][1]
   The condition that ('a', 'b') == (1, 3), (3, 1), (5, 7), (7, 5) mod 8

We also show that the result presented in proposition 3.6 is correct.

::

   sage: Ec2 = Ec.twist(2)
   sage: C2 = C & CongruenceCondition(a + b, 2) & ~CongruenceCondition(a + b, 4)
   sage: Ec2.conductor_exponent(B2, condition=C2)
   4 if ('a', 'b') is 1 of 4 possibilities mod 8
   0 if ('a', 'b') is 1 of 4 possibilities mod 8

We compute the conductor of the restriction of scalars as is done in
proposion 4.1 and proposition 4.2.

::

   sage: Pbad = Ec.decomposition_field().primes_above(2*5)
   sage: Ec.conductor_restriction_of_scalars(additive_primes=Pbad)
   2^(2*n0+8)*5^(n1+6)*Norm(Rad_P( ((49280*zeta0^3 - 130240*zeta0^2 - 41600*zeta0 + 211200)) * (a^2 + (-zeta0^2 + 2)*a*b + b^2) * (a^2 + (zeta0^2 - 3)*a*b + b^2)^2 ))
    where 
   n0 =  8 if ('a', 'b') is 1 of 6 possibilities mod 4
         6 if ('a', 'b') is 1 of 4 possibilities mod 4
	 4 if ('a', 'b') is 1 of 4 possibilities mod 8
	 0 if ('a', 'b') is 1 of 4 possibilities mod 8
   n1 =  2 if ('a', 'b') is 1 of 20 possibilities mod 5
         0 if ('a', 'b') is 1 of 4 possibilities mod 5

The levels of the newforms is as mentioned at the end of section 4 one
of 100, 4000, 800 or 1600. Note that our computation also allows other
(lower) levels.

::

   sage: Ec.newform_levels(bad_primes=Pbad)
   [(1600, 1600)]             if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(320, 1600), (1600, 320)] if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(800, 800)]               if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(160, 800), (800, 160)]   if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(400, 400)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(80, 400), (400, 80)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(100, 100)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(20, 100), (100, 20)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5

We circumvent the code choosing the lower levels instead of the levels
we want by explicitly computing the spaces of newforms as in the
article.

::

   sage: level = apply_to_conditional_value(lambda ls: ls[0][1], Ec.newform_levels(bad_primes=Pbad))
   sage: char = Ec.splitting_character('conjugacy')[1]^(-1)
   sage: nfs = apply_to_conditional_value(lambda lvl: get_newforms(lvl, character=char,
   ....: algorithm='sage'), level)

As in the article we divide these spaces into three different
categories ``S1``, ``S2`` and ``S3``, respectively the newforms with
complex multiplication, those without CM and a coefficient field of
degree strictly larger than 2, and those without CM and coefficient
field :math:`\QQ(i)`.

::

  sage: S1 = apply_to_conditional_value(lambda ls: [nf for nf in ls if nf.has_cm()], nfs)
  sage: S2 = apply_to_conditional_value(lambda ls: [nf for nf in ls if not nf.has_cm() and
  ....: nf.coefficient_field().absolute_degree() > 2], nfs)
  sage: S3 = apply_to_conditional_value(lambda ls: [nf for nf in ls if not nf.has_cm() and
  ....: nf.coefficient_field().is_isomorphic(QQ[sqrt(-1)])], nfs)

Case 2 divides d
----------------
  
The article reasons that in this case the level 800 does not
appear. The article claims there are 8 newforms in ``S1`` of which
half have complex multiplication by :math:`\QQ(i)` and the other half
have complex multiplication by :math:`\QQ(\sqrt{5})`.

::

   sage: len(S1[0][0] + S1[2][0] + S1[3][0])
   10
   sage: [nf._f.cm_discriminant() for nf in S1[0][0] + S1[2][0] + S1[3][0]]
   [-20, -4, -4, -20, -20, -4, -20, -20, -20, -4]

These newforms are eliminated for all primes :math:`p > 13` for which
:math:`p \equiv 1` modulo 4 or :math:`p \equiv \pm 1` modulo 5.

The article now claims there are 12 newforms in ``S2``.

::

   sage: len(S2[0][0] + S2[1][0] + S2[2][0])
   13

Furthermore it eliminates all these newforms by the fact that none of
them have a third coefficient of the form :math:`t - i t` with
:math:`t` an integer modulo some prime ``P`` above :math:`p`, whenever
:math:`p > 5`. We compute all prime numbers below primes that could
divide the difference between the third coefficient of a newform in
``S2`` and :math:`t - i t` for all possible :math:`|t| \le 2` as
mentioned in the article.

::

   sage: lcm(ZZ((nf.coefficient(3) - t*(1 - sqrt(nf.coefficient_field()(-1)))).absolute_norm())
   ....: for nf in S2[0][0] + S2[2][0] + S2[3][0] for t in range(-2, 3)).prime_factors()
   [2, 3, 5, 7, 29]

As claimed in the article we check there is 10 newforms in ``S3``.

::

   sage: len(S3[0][0])
   10

The article claims that the twist of each newform in ``S3`` by the
character of :math:`\QQ(sqrt{2})` is a newform of level 800, which we
check.

::

   sage: chi = character_for_root(2)
   sage: all(any(all(nf.coefficient(i) * chi(i) == ng.coefficient(i) for i in range(sturm_bound(800)+1))
   ....: for ng in S3[1][0]) for nf in S3[0][0])
   True

Next the article has an argument to prove that this is impossible for
newforms associated to our problem eliminating all newforms from
``S3``.

Case 3 divides d
----------------

We first check the claim the article makes about the additional
newforms in this case at level 800, namely 0 in ``S1``, 4 in ``S2``
and 10 in ``S3``.

::

   sage: len(S1[1][0])
   0
   sage: len(S2[1][0])
   4
   sage: len(S3[1][0])
   10

The article looks first at the case when :math:`a + b` is odd, in
which case only the levels 800 and 1600 are relevant. For the newforms
of level 800 in ``S2`` the article appplies the same trick as
before. In this case they note that :math:`t - i t` for :math:`t` an
integer of absolute value at most 2 is not a root of the minimal
polynomial of the third coefficient of a newform over :math:`\QQ(i)`
modulo primes above :math:`p > 73`. We compute all the prime numbers
below primes dividing :math:`t - i t` substituted in such a minimal
polynomial to verify this.

::

   sage: i = nf.coefficient_field().base_field().gen()
   sage: lcm(ZZ(nf.coefficient(3).minpoly()(t - t*i).norm()) for nf in S2[1][0]
   ....: for t in range(-2, 3)).prime_factors()
   [5, 29]

The newforms of level 800 in ``S3`` are eliminated by comparing traces
of frobenius which we verify. We verify that only for primes above 2,
3 and 5 these traces can be the same.

::

   sage: result = eliminate_by_trace(Ec, S3[1][0], 3, condition=C & CongruenceCondition(a + b, 3))
   sage: lcm(nf[1] for nf in result).prime_factors()
   [2, 3, 5]

Next, the article turns to the newforms of level 1600 in the case when
:math:`a + b` is odd. For those in ``S1`` the article notes that those
with complex multiplication by :math:`\QQ(\sqrt{-5})` can be
eliminated in the way we eliminated newforms in ``S3`` of
level 800. The others can be eliminated using the same argument as for
:math:`2 \mid d`. We first remove the latter and show that the primes
for which the comparison of traces at 3 could still work are at
most 5.

::

   sage: result = [nf for nf in S1[0][0] if nf._f.cm_discriminant() == -20]
   sage: result = eliminate_by_trace(Ec, result, 3, condition=C & CongruenceCondition(a + b, 3))
   sage: lcm(nf[1] for nf in result).prime_factors()
   [2, 3, 5]

The article reasons that for the newforms of level 1600 in ``S2`` the
same argument as in the case :math:`2 \mid d` holds. For the remaining
newforms of level 1600 in ``S3`` the article uses a similar argument
as in the case :math:`2 \mid d`. This concludes all cases with
:math:`a + b` odd.

Now for the case that :math:`a + b` is even the article reasons the
previously computed results are sufficient to reduce to the cases for
theorem 5.2.

Multi-Frey approach
-------------------

For the full results the article uses a new Frey curve which is also a
Q-curve.

::

   sage: a_invariants2 = [0, 2*(a - b), 0, (-3/10*sqrt(K(5)) + 1/2)*phi1, 0]
   sage: isogenies2 = {sigma^0: (QQ(1), 1), sigma^1: (sqrtm2, 2)}
   sage: F = FreyQcurve(a_invariants2, isogenies=isogenies2, condition=C)

The article claims that :math:`F` has the same splitting behaviour as
:math:`E` and that twisting by the same :math:`\gamma` gives a
decomposable twist, which we check.

::

   sage: E.splitting_character() == F.splitting_character()
   True
   sage: E.splitting_field().is_isomorphic(F.splitting_field())
   True
   sage: F.does_decompose()
   False
   sage: Fc = F.twist(gamma)
   sage: Fc.does_decompose()
   True

According to the article the corresponding newforms have level 100 if
:math:`8 \mid a + b`, 400 if :math:`4 \| a + b` or 1600 if :math:`2 \|
a + b`, which we check.

::

   sage: Fc.newform_levels()
   Warning: Assuming that a and b are coprime.
   [(320, 1600), (1600, 320)] if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(1600, 1600)]             if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(160, 800), (800, 160)]   if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(800, 800)]               if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(80, 400), (400, 80)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(400, 400)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(20, 100), (100, 20)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(100, 100)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5

Our code produces some lower levels than the levels mentioned in the
article, hence to stick to the levels mentioned in the article we
omit the newform_candidates method.

::

   sage: levels2 = apply_to_conditional_value(lambda ls: ls[0][1], Fc.newform_levels())
   Warning: Assuming that a and b are coprime.
   sage: char2 = Fc.splitting_character('conjugacy')[1]^(-1)
   sage: nfs2 = apply_to_conditional_value(lambda lvl: get_newforms(lvl, character=char,
   ....: algorithm='sage'), levels2)

The article remarks that all the pairs :math:`(f, g)` of newforms, one
for ``Ec`` and one for ``Fc`` respectively, for which :math:`f` does
not have CM can be removed by previous arguments. Similarly can those
for which :math:`f` or :math:`g` has a coefficient field strictly
larger than :math:`\QQ(i)`. As in the article we apply multi-Frey
comparison of traces at 3, 7, 13 and 17.

::

   sage: nfs22 = apply_to_conditional_value(lambda ls: [nf for nf in ls
   ....: if nf.coefficient_field().absolute_degree() == 2], nfs2)
   sage: S12 = apply_to_conditional_value(lambda ls: [nf for nf in ls
   ....: if nf.coefficient_field().absolute_degree() == 2], S1)
   sage: nfs_big = conditional_product(S1, nfs22)
   sage: nfs_big = ConditionalValue([(val, con) for val, con in nfs_big
   ....: if not con.pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()])
   sage: nfs_big = eliminate_by_traces((Ec, Fc), nfs_big, primes=[3, 7, 13, 17])

According to the article we should only have 8 newforms remaining if
we remove all cases in which only a prime :math:`p \le 13` would work.

::

   sage: nfs_big = eliminate_primes((Ec, Fc), nfs_big, product(prime_range(14)))
   sage: sum(len(nfs_big[i][0]) for i in range(len(nfs_big)))

We however find there is more newforms remaining than only CM forms,
but removing the case :math:`p = 17` we are indeed in the case as
described that all pairs :math:`(f, g)` have CM with different
discriminants.

::

   sage: nfs_big = eliminate_primes((Ec, Fc), nfs_big, 17)
   sage: sum(len(nfs_big[i][0]) for i in range(len(nfs_big)))
   12
   sage: [fg[0]._f.cm_discriminant() != fg[1]._f.cm_discriminant()
   ....:  for fg in nfs_big[0][0] + nfs_big[1][0] + nfs_big[2][0] + nfs_big[3][0]]
   [True, True, True, True, True, True, True, True, True, True, True, True]

There are however more newforms as expressed in the article.
