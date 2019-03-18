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
   sage: (2*w + 1)^2
   5
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

   sage: Ec.conductor_restriction_of_scalars(additive_primes=Ec.decomposition_field().primes_above(2*5))
   2^(2*n0+8)*5^(n1+6)*Norm(Rad_P( ((49280*zeta0^3 - 130240*zeta0^2 - 41600*zeta0 + 211200)) * (a^2 + (-zeta0^2 + 2)*a*b + b^2) * (a^2 + (zeta0^2 - 3)*a*b + b^2)^2 ))
    where 
   n0 =  8 if ('a', 'b') is 1 of 6 possibilities mod 4
         6 if ('a', 'b') is 1 of 4 possibilities mod 4
	 4 if ('a', 'b') is 1 of 4 possibilities mod 8
	 0 if ('a', 'b') is 1 of 4 possibilities mod 8
   n1 =  2 if ('a', 'b') is 1 of 20 possibilities mod 5
         0 if ('a', 'b') is 1 of 4 possibilities mod 5

The levels of the newforms is as mentioned at the end of section 5 one
of 100, 4000, 800 or 1600. Note that our computation also allows other
(lower) levels.

::

   sage: Ec.newform_levels(bad_primes=Ec.decomposition_field().primes_above(2*5))
   [(1600, 1600)]             if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(320, 1600), (1600, 320)] if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(800, 800)]               if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(160, 800), (800, 160)]   if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(400, 400)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(80, 400), (400, 80)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5
   [(100, 100)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
   [(20, 100), (100, 20)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5

