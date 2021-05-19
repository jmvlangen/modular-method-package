=======================================
 Second returning example of Chapter 2
=======================================

This file contains the code for the second returning example of
Chapter 2.

.. linkall

Example 2.9.10
--------------

We start by defining the curve and check that the decomposition field
is as claimed in the example.

::

   sage: L.<zeta40> = CyclotomicField(40)
   sage: K.<t> = L.subfield(zeta40 + zeta40^(-1))[0]
   sage: sqrt2, sqrt5 = sqrt(K(2)), sqrt(K(5))
   sage: c = sqrt((5 + sqrt5) / 2)
   sage: gamma = (1 + sqrt2)*(sqrt5 + c)
   sage: E = Qcurve([0, 4*gamma, 0, 2*gamma^2*(1 + sqrt2/sqrt5), 0],
   ....:            guessed_degrees=[2])
   sage: E.decomposition_field() == K
   True
   sage: E.does_decompose()
   True

We compute a twist character and a splitting image field for each
Galois orbit of splitting maps.

::

   sage: E.twist_character('conjugacy')
   (Dirichlet character modulo 1 of conductor 1,
    Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4)
   sage: E.splitting_image_field('conjugacy')
   (Cyclotomic Field of order 8 and degree 4,
    Cyclotomic Field of order 8 and degree 4)

Now we compute the right hand side of Equation (2.7).

::

   sage: RHS = E.conductor_restriction_of_scalars()
   sage: RHS.factor()
   2^72 * 3^8 * 5^12

We compute the splitting characters as in the example.

::

   sage: E.splitting_character('conjugacy')
   (Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4,
    Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -zeta4)

Finally we verify the newform levels computed in the example.

::

   sage: E.newform_levels()
   [(7680, 38400), (38400, 7680)]


Example 2.10.8
--------------

We verify that the curve has multiplicative reduction at a prime above 3.

::

   sage: P3 = K.prime_above(3)
   sage: E.has_multiplicative_reduction(P3)
   True

We compute :math:`-c_4 / c_6` and check it is not a square in `K`.

::

   sage: gamma = -E.c4() / E.c6()
   sage: gamma.is_square()
   False

We show that 3 does not ramify in the extension
:math:`K(\sqrt{\gamma})`.

::

   sage: R.<x> = QQ[]
   sage: L.<a> = gamma.minpoly()(x^2).splitting_field()
   sage: G = L.galois_group()
   sage: iota = K.embeddings(L)[0]
   sage: sqrtgamma = sqrt(iota(gamma))
   sage: Q3 = Kext.prime_above(iota(P3))
   sage: Q3.ramification_index()
   1

We compute the Artin symbol of this prime above 3 and check it is the
same for all primes above 3.

::

   sage: s3 = G.artin_symbol(Q3)
   sage: all(s3 == G.artin_symbol(Q) for Q in L.primes_above(3))
   True

Now we perform the computation described in Theorem 2.10.7 and find
the results as mentioned in the article.

::

   sage: from modular_method.number_fields.galois_group import galois_field_restrict
   sage: s3K = galois_field_restrict(s3, K, embedding=iota)
   sage: as3 = QQ(iota(E.isogeny_scalar(s3K)) * sqrtgamma / s3(sqrtgamma)); as3
   2
   sage: Lbeta = E.splitting_image_field(); zeta8 = Lbeta.gen()
   sage: trace = E.splitting_map()(s3K)^(-1) * as3
   sage: trace == zeta8^3
   True
