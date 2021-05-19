=================================================
 Returning example from Chapter 2 (Frey variant)
=================================================

This file contains the code for the Frey Q-curve of the returning
example from Chapter 2.

.. linkall

Example 2.1.4
-------------

We enter the Frey Q-curve in the framework.

::

   sage: R.<a,b> = QQ[]
   sage: _.<sqrt3> = QuadraticField(3)
   sage: con = CoprimeCondition([a, b])
   sage: E = FreyQcurve([0, 12*a, 0, 18*(a^2 + b*sqrt3), 0],
   ....:                condition=con, guessed_degrees=[2]); E
   Frey Q-curve defined by y^2 = x^3 + 12*a*x^2 + (18*a^2+(18*sqrt3)*b)*x over Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878? with parameters (a, b)

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
   sage: [E.isogeny_scalar(s)^2 for s in [G(1), s2, s3, s2*s3]]
   [1, -2, 1, -2]

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
   Frey Q-curve defined by y^2 = x^3 + ((-6*lu0-12)*a)*x^2 + ((18*lu0+72)*a^2+(36*lu0+108)*b)*x over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu with parameters (a, b)

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

   sage: RHS = Egamma.conductor_restriction_of_scalars(); RHS
   2^(n0+4)*3^(n1+2)*Norm(Rad_P( ((-22394880*lu0 + 77635584)) * (a^2 + (-1/2*lu0)*b) * (a^2 + (1/2*lu0)*b)^2 ))
    where
   n0 =  12 if ('a', 'b') == (1, 0) mod 2
         14 if ('a', 'b') == (1, 1) mod 2
         8  if ('a', 'b') == (0, 3), (2, 3) mod 4
         0  if ('a', 'b') is 1 of 4 possibilities mod 8
         4  if ('a', 'b') is 1 of 4 possibilities mod 8
   n1 =  0 if ('a', 'b') is 1 of 6 possibilities mod 3
         2 if ('a', 'b') == (0, 1), (0, 2) mod 3

Example 2.10.4
--------------

Still to be worked out. The code below is a copy from the code for a
regular Q-curve.

::

   sage: sqrt3 = K.gen()
   sage: Egamma = E.twist(gamma)
   sage: Egood = Egamma.global_minimal_model()
   sage: Kgood = Egood.base_ring()
   sage: Egood = Egood.change_ring(Kgood.hom([-2*sqrt3]))
   sage: Egood = Qcurve(Egood.rst_transform(-1 - sqrt3, 0, 0), guessed_degrees=[2]); Egood
   Q-curve defined by y^2 = x^3 + (-2*sqrt3-2)*x^2 + (3*sqrt3+5)*x over Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?
   sage: [P.smallest_integer() for P, e in (K.ideal(Egood.c4()) + K.ideal(Egood.discriminant())).factor()]
   [2]
   sage: [Egood.isogeny_scalar(s) for s in K.galois_group()]
   [1, -sqrt3 - 1]
   sage: [Egood.isogeny_x_map(s) for s in K.galois_group()]
   [x, ((-1/2*sqrt3 + 1)*x^2 + (-sqrt3 + 1)*x + 1/2*sqrt3 + 1/2)/x]
   sage: R.<x> = K[]
   sage: f = x^2 + Egood.a2()*x + Egood.a4()
   sage: F = ((2 - sqrt3)/2) * (f / x)
   sage: Egood.isogeny_x_map(K.galois_group().gen()) == F
   True
   sage: R = 4*x^3 + Egood.b2()*x^2 + 2*Egood.b4()*x + Egood.b6()
   sage: R == 4*x*f
   True
   sage: Lbeta.<sqrtm2> = QuadraticField(-2)
   sage: def mylatex(n):
   ....:     return "${ " + latex(n) + " }$"
   sage: for p in prime_range(3, 30):
   ....:     P = K.prime_above(p)
   ....:     FP = P.residue_field()
   ....:     if FP.degree() == 1:
   ....:         # The case sigma in G_K
   ....:         trace = 1 + p - Egood.reduction(P).cardinality()
   ....:         print(mylatex(p), "&", "ramifies/splits", "&", mylatex(trace), "\\\\")
   ....:     if len(K.primes_above(p)) == 1:
   ....:         # The case sigma not in G_K
   ....:         c1_ = 2*x^(p + 1) - (2 - sqrt3)*f
   ....:         c2_ = 2^p * (1 - sqrt3) * x^((p + 3)/2) * f^((p - 1)/2) - x^2 + 5 + 3*sqrt3
   ....:         c1_ = c1_.change_ring(FP)
   ....:         c2_ = c2_.change_ring(FP)
   ....:         trace = sqrtm2 * (gcd(c1_, c2_).radical().degree() - (p + 1)/2)
   ....:         print(mylatex(p), "&", "ramifies/inert", "&", mylatex(trace), "\\\\")

Example 2.10.10
---------------

We verify the results obtained in this example.

::

   sage: Egamma.trace_of_frobenius(7)
   0          if ('a', 'b') is 1 of 12 possibilities mod 7
   -3*zeta4a0 if ('a', 'b') is 1 of 6 possibilities mod 7
   zeta4a0    if ('a', 'b') is 1 of 6 possibilities mod 7
   -2*zeta4a0 if ('a', 'b') is 1 of 6 possibilities mod 7
   2*zeta4a0  if ('a', 'b') is 1 of 6 possibilities mod 7
   -zeta4a0   if ('a', 'b') is 1 of 6 possibilities mod 7
   3*zeta4a0  if ('a', 'b') is 1 of 6 possibilities mod 7
   sage: Egamma.trace_of_frobenius(7, splitting_map=1)
   0          if ('a', 'b') is 1 of 12 possibilities mod 7
   3*zeta4a0  if ('a', 'b') is 1 of 6 possibilities mod 7
   -zeta4a0   if ('a', 'b') is 1 of 6 possibilities mod 7
   2*zeta4a0  if ('a', 'b') is 1 of 6 possibilities mod 7
   -2*zeta4a0 if ('a', 'b') is 1 of 6 possibilities mod 7
   zeta4a0    if ('a', 'b') is 1 of 6 possibilities mod 7
   -3*zeta4a0 if ('a', 'b') is 1 of 6 possibilities mod 7
   sage: Egamma.trace_of_frobenius(11)
   6   if ('a', 'b') is 1 of 5 possibilities mod 11
   -6  if ('a', 'b') is 1 of 5 possibilities mod 11
   0   if ('a', 'b') is 1 of 30 possibilities mod 11
   -4  if ('a', 'b') is 1 of 20 possibilities mod 11
   -2  if ('a', 'b') is 1 of 10 possibilities mod 11
   2   if ('a', 'b') is 1 of 10 possibilities mod 11
   4   if ('a', 'b') is 1 of 20 possibilities mod 11
   12  if ('a', 'b') is 1 of 10 possibilities mod 11 and a11E == +1 or ('a', 'b') is 1 of 10 possibilities mod 11 and a11E == +1
   -12 if ('a', 'b') is 1 of 10 possibilities mod 11 and a11E == -1 or ('a', 'b') is 1 of 10 possibilities mod 11 and a11E == -1
