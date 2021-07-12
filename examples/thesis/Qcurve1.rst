===================================
 Returning example from Chapter 2
===================================

This file contains the code for the Q-curve of the returning example
from Chapter 2.

.. linkall

The following import is required for all examples to work

::

   sage: from modular_method import *

Example 2.1.4
-------------

We enter the Q-curve in the framework.

::

   sage: _.<sqrt3> = QuadraticField(3)
   sage: E = Qcurve([0, 12, 0, 18*(1 + sqrt3), 0],
   ....:            guessed_degrees=[2]); E
   Q-curve defined by y^2 = x^3 + 12*x^2 + (18*sqrt3+18)*x over Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?

Example 2.2.3
-------------

We verify the results of this example. First we do the isogeny
computation as in the article.

::

   sage: K = E.complete_definition_field()
   sage: K.is_isomorphic(QQ[sqrt(-2), sqrt(3)])
   True
   sage: sqrtm2, sqrt3 = sqrt(K(-2)), sqrt(K(3))
   sage: G = K.galois_group()
   sage: s2 = next(s for s in G if s != G(1) and s(sqrtm2) == sqrtm2)
   sage: s3 = next(s for s in G if s != G(1) and s(sqrt3) == sqrt3)
   sage: (E.isogeny_scalar(G(1)) == 1 and
   ....:  E.isogeny_scalar(s2) == sqrtm2 and
   ....:  E.isogeny_scalar(s3) == 1 and
   ....:  E.isogeny_scalar(s2*s3) == sqrtm2)
   True

Next we verify that the scalars of the corresponding isogenies are
indeed as mentioned in the example.

::

   sage: test = {s: E.degree_map(s) / E.isogeny_scalar(s) for s in [G(1), s2, s3, s2*s3]}
   sage: (test[G(1)] == 1 and
   ....:  test[s2] == -sqrtm2 and
   ....:  test[s3] == 1 and
   ....:  test[s2*s3] == -sqrtm2)
   True

Lastly we check that the table for the coboundary :math:`c_E` is
correct.

::

   sage: matrix([[E.c(s, t)
   ....:          for t in [G(1), s2, s3, s2*s3]]
   ....:         for s in [G(1), s2, s3, s2*s3]])
   [ 1  1  1  1]
   [ 1 -2  1 -2]
   [ 1 -1  1 -1]
   [ 1  2  1  2]

Example 2.4.4
-------------

We verify the results mentioned in this example.

::

   sage: [E.degree_map(s) for s in [G(1), s2, s3, s2*s3]]
   [1, 2, 1, 2]
   sage: E.dual_basis()
   ([3], [2])
   sage: E.xi_pm()
   [(3, 2)]
   sage: E.xi_pm_local(2), E.xi_pm_local(3), E.xi_pm_local(5)
   (-1, -1, 1)
   sage: eps = E.splitting_character()
   sage: eps == next(eps for eps in DirichletGroup(12)
   ....:             if eps.conductor() == 12)
   True
   sage: beta = E.splitting_map()
   Warning: The restriction of scalars of this Q-curve over the decomposition field does not decompose into abelian varieties of GL_2-type. Use the method decomposable_twist to find a twist that does.
   sage: [beta(s)^2 for s in [G(1), s2, s3, s2*s3]]
   [1, -2, 1, -2]

Example 2.5.4
-------------

First of all we verify that the matrix given in this example is
correct.

::

   sage: matrix([[E.c_splitting_map(s, t) / E.c(s, t)
   ....:          for t in [G(1), s2, s3, s2*s3]]
   ....:         for s in [G(1), s2, s3, s2*s3]])
   [ 1  1  1  1]
   [ 1  1  1  1]
   [ 1 -1  1 -1]
   [ 1 -1  1 -1]

Next we show that the map :math:`\alpha` given in the example indeed
has coboundary :math:`c_\beta c_E^{-1}`.

::

   sage: alpha = {s : 1 if s in [G(1), s3] else (1 + sqrt3)/sqrtm2
   ....:          for s in [G(1), s2, s3, s2*s3]}
   sage: all(E.c_splitting_map(s, t) / E.c(s, t) == alpha[s]*s(alpha[t])*(alpha[s*t])^(-1)
   ....:     for s in G for t in G)
   True

We check that the :math:`\gamma` computed from :math:`\alpha` indeed
satisfies the claimed properties.

::

   sage: gamma = sum(alpha[s]^(-2) for s in [G(1), s2, s3, s2*s3])
   sage: gamma == -2 + 2*sqrt3
   True
   sage: all(alpha[s]^2 == s(gamma) / gamma for s in [G(1), s2, s3, s2*s3])
   True

We verify the computation of the table for :math:`c_\beta c_E^{-1}`
over :math:`K_\gamma`.

::

   sage: gamma = 1 - sqrt3
   sage: R.<x> = K[]
   sage: Kgamma.<sqrtgamma> = K.extension(x^2 - gamma)
   sage: sqrtm6 = Kgamma(sqrtm2*sqrt3)
   sage: Kgamma.<a> = Kgamma.absolute_field()
   sage: sqrtgamma, sqrtm6 = Kgamma(sqrtgamma), Kgamma(sqrtm6)
   sage: Ggamma = Kgamma.galois_group()
   sage: sgamma = next(s for s in Ggamma
   ....:               if s != Ggamma(1) and
   ....:               s(sqrtgamma) == sqrtgamma)
   sage: s6 = next(s for s in Ggamma
   ....:           if s(sqrt(Kgamma(-2))) != sqrt(Kgamma(-2)) and
   ....:           s(sqrtm6) == sqrtm6)
   sage: Gls = [Ggamma(1), s6, s6^2, s6^3,
   ....:        sgamma, s6*sgamma, s6^2*sgamma, s6^3*sgamma]
   sage: all(s in Gls for s in Ggamma)
   True
   sage: matrix([[E.c_splitting_map(s, t) / E.c(s, t) for t in Gls] for s in Gls])
   [ 1  1  1  1  1  1  1  1]
   [ 1 -1  1 -1  1 -1  1 -1]
   [ 1  1  1  1  1  1  1  1]
   [ 1 -1  1 -1  1 -1  1 -1]
   [ 1 -1  1 -1  1 -1  1 -1]
   [ 1  1  1  1  1  1  1  1]
   [ 1 -1  1 -1  1 -1  1 -1]
   [ 1  1  1  1  1  1  1  1]

Next we do the confirmation of the map :math:`\alpha`.

::

   sage: alpha = {s : 1 if s in [Ggamma(1), s6, sgamma, s6*sgamma] else -1
   ....:          for s in Gls}
   sage: all(E.c_splitting_map(s, t) / E.c(s, t) ==
   ....:     alpha[s] * alpha[t] / alpha[s*t]
   ....:     for s in Gls for t in Gls)
   True

We check the splitting map as in the example, but also confirm this is
the same as the one computed here.

::

   sage: beta = {s : E.splitting_map()(s) * alpha[s]
   ....:         for s in Gls}
   sage: betasqrtm2 = E.splitting_image_field().gen()
   sage: (betasqrtm2^2 == -2 and
   ....:  beta[Ggamma(1)] == 1 and
   ....:  beta[s6] == betasqrtm2 and
   ....:  beta[s6^2] == -1 and
   ....:  beta[s6^3] == -betasqrtm2 and
   ....:  beta[sgamma] == 1 and
   ....:  beta[s6*sgamma] == betasqrtm2 and
   ....:  beta[s6^2*sgamma] == -1 and
   ....:  beta[s6^3*sgamma] == -betasqrtm2)
   True
   sage: all(E.c(s, t) == beta[s] * beta[t] / beta[s*t] for s in Gls for t in Gls)
   True

Example 2.6.1
-------------

We confirm that there are four splitting maps and the corresponding
non-trivial twist characters are the quadratic characters of
:math:`\QQ(\sqrt{-2})`, :math:`\QQ(\sqrt{3})`, and
:math:`\QQ(\sqrt{-6})`.

::

   sage: iota = E.definition_field().embeddings(Kgamma)[0]
   sage: Egamma = E.change_ring(iota)
   sage: Egamma.number_of_splitting_maps()
   4
   sage: chis = Egamma.twist_character('all', galois=True)
   sage: kernels = [Ggamma.subgroup(s for s in Ggamma if chi(s) == 1)
   ....:            for chi in chis]
   sage: fields = [kernel.fixed_field()[0] for kernel in kernels]
   sage: [(field.degree(), field.discriminant().squarefree_part())
   ....:  for field in fields]
   [(1, 1), (2, -2), (2, 3), (2, -6)]

Next we compute the number of splitting maps and one splitting map
within each Galois orbit.

::

   sage: Egamma.number_of_splitting_maps(count_conjugates=False)
   2
   sage: beta1, beta2 = Egamma.splitting_map('conjugacy')
   sage: Lbeta = Egamma.splitting_image_field()
   sage: Gbeta = Lbeta.galois_group()
   sage: all(any(beta1(s) != t(beta2(s)) for s in Ggamma)
   ....:     for t in Gbeta)
   True

We verify that these splitting maps agree with the ones given in the
example.

::

   sage: (beta1(Ggamma(1)) == 1 and
   ....:  beta1(s6) == -betasqrtm2 and
   ....:  beta1(s6^2) == -1 and
   ....:  beta1(s6^3) == betasqrtm2 and
   ....:  beta1(sgamma) == 1 and
   ....:  beta1(s6*sgamma) == -betasqrtm2 and
   ....:  beta1(s6^2*sgamma) == -1 and
   ....:  beta1(s6^3*sgamma) == betasqrtm2 and
   ....:  beta2(Ggamma(1)) == 1 and
   ....:  beta2(s6) == betasqrtm2 and
   ....:  beta2(s6^2) == -1 and
   ....:  beta2(s6^3) == -betasqrtm2 and
   ....:  beta2(sgamma) == -1 and
   ....:  beta2(s6*sgamma) == -betasqrtm2 and
   ....:  beta2(s6^2*sgamma) == 1 and
   ....:  beta2(s6^3*sgamma) == betasqrtm2)
   True

Example 2.7.9
-------------

We compute the degree field.

::

   sage: E.degree_field()
   Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?

Next we compute the twist of the curve discussed in the example.

::

   sage: E.decomposable_twist()
   Q-curve defined by y^2 = x^3 + (-6*lu0-12)*x^2 + (-18*lu0-36)*x over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu

Example 2.9.3
-------------

First of all we perform the twist on the curve.

::

   sage: Egamma = E.twist(gamma)

Next we verify that the splitting image field is indeed
:math:`\Q(\sqrt{-2})`

::

   sage: Egamma.splitting_image_field().is_isomorphic(QuadraticField(-2))
   True

Now we compute the conductor of the restriction of scalars for
`Egamma`.

::

   sage: RHS = Egamma.conductor_restriction_of_scalars()
   sage: RHS.factor()
   2^18 * 3^2

Example 2.10.4
--------------

We start by computing a global minimal model of :math:`E_\gamma` and
verify it is the same as the one given in the example.

::

   sage: Kgood = Egamma.definition_field()
   sage: sqrt3 = sqrt(Kgood(3))
   sage: Egood = Qcurve(Egamma.scale_curve(1/2 + sqrt3/6), guessed_degrees=[2])
   sage: Egood.is_global_minimal_model()
   True
   sage: (Egood.a1() == 0 and
   ....:  Egood.a2() == -2*(1 + sqrt3) and
   ....:  Egood.a3() == 0 and
   ....:  Egood.a4() == -1*(1 + sqrt3) and
   ....:  Egood.a6() == 0)
   True

Next we show that :math:`c_4` and the discriminant of this curve are
coprime outside primes above 2.

::

   sage: [P.smallest_integer() for P, e in (K.ideal(Egood.c4()) + K.ideal(Egood.discriminant())).factor()]
   [2]

We verify the invariants of the isogenies in the example are correct.

::

   sage: Ggood.<s3> = Kgood.galois_group()
   sage: _.<x> = Kgood[]
   sage: f = x^2 + Egood.a2()*x + Egood.a4()
   sage: F = ((2 - sqrt3)/2) * (f / x)
   sage: (Egood.isogeny_scalar(Ggood(1)) == 1 and
   ....:  Egood.isogeny_x_map(Ggood(1)) == x and
   ....:  Egood.isogeny_scalar(s3) == -1 - sqrt3 and
   ....:  Egood.isogeny_x_map(s3) == F)
   True

We compute the polynomial :math:`R` and verify :math:`f_1` and
:math:`f_2` are correct by computing the numerator and denominator of
both :math:`F(x)` and :math:`F'(x) R`.

::

   sage: R = 4*x^3 + Egood.b2()*x^2 + 2*Egood.b4()*x + Egood.b6()
   sage: (R == 4*x*f and
   ....:  F.numerator() == ((2 - sqrt3) / 2) * f and
   ....:  F.denominator() == x and
   ....:  (F.derivative() * R).numerator() == 2 * (2 - sqrt3) * f * (x^2 + 1 + sqrt3) and
   ....:  (F.derivative() * R).denominator() == x)
   True

Finally we verify all the values in Table 2.1.

::

   sage: Lbeta.<sqrtm2> = QuadraticField(-2)
   sage: for p in prime_range(3, 30):
   ....:     P = Kgood.prime_above(p)
   ....:     FP = P.residue_field()
   ....:     if FP.degree() == 1:
   ....:         # The case sigma in G_K
   ....:         trace = 1 + p - Egood.reduction(P).cardinality()
   ....:         print(p, "ramifies/splits, trace:", trace)
   ....:     if len(Kgood.primes_above(p)) == 1:
   ....:         # The case sigma not in G_K
   ....:         f1_ = 2*x^(p + 1) - (2 - sqrt3)*f
   ....:         f2_ = 2^p * (1 + sqrt3) * x^((p + 3)/2) * f^((p - 1)/2) + (2 - sqrt3) * (x^2 + 1 + sqrt3)
   ....:         f1_ = f1_.change_ring(FP)
   ....:         f2_ = f2_.change_ring(FP)
   ....:         trace = sqrtm2 * (gcd(f1_, f2_).radical().degree() - (p + 1)/2)
   ....:         print(p, "ramifies/inert, trace:", trace)
   3 ramifies/splits, trace: -2
   3 ramifies/inert, trace: 2*sqrtm2
   5 ramifies/inert, trace: sqrtm2
   7 ramifies/inert, trace: -3*sqrtm2
   11 ramifies/splits, trace: -4
   13 ramifies/splits, trace: -2
   17 ramifies/inert, trace: 2*sqrtm2
   19 ramifies/inert, trace: 0
   23 ramifies/splits, trace: 8
   29 ramifies/inert, trace: 5*sqrtm2
