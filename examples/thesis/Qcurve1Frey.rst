=================================================
 Returning example from Chapter 2 (Frey variant)
=================================================

This file contains the code for the Frey Q-curve of the returning
example from Chapter 2.

.. linkall

The following import is required for all examples to work

::

   sage: from modular_method import *

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
   sage: s2 = next(s for s in G if s != G(1) and
   ....:           s(sqrtm2) == sqrtm2)
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
   Warning: Assuming that (-1440*lu0 + 5760)*a^2 + (-1728*lu0 + 5184)*b and (-22394880*lu0 + 77635584)*a^6 + (38817792*lu0 - 134369280)*a^4*b + (67184640*lu0 - 232906752)*a^2*b^2 + (-116453376*lu0 + 403107840)*b^3 are coprime outside ('(1/2*lu0 + 1)', '(1/2*lu0)').
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

First of all we determine the primes for which the invariants
:math:`c_4` and :math:`Delta` of `Egamma` are not coprime by computing
their resultant.

::

   sage: Kgood = Egamma.definition_field()
   sage: sqrt3 = sqrt(Kgood(3))
   sage: set(P.smallest_integer()
   ....:     for P, _ in Kgood.ideal(Egamma.c4().resultant(Egamma.discriminant())(1, 1)).factor())
   {2, 3}

We thus see that the model of :math:`Egamma` is minimal at all primes
not above 2 and 3. Furthermore the :math:`Egamma` has good reduction
at such a prime `P` if `P` does not divide the discriminant, and
multiplicative reduction otherwise. We calculate the possible
reductions types at primes above 2 and 3 by computing the conductor
exponents at those primes. Note that 2 and 3 ramify in `Kgood` so
there is only one prime above each of them.

::
   
   sage: Kgood.discriminant().prime_factors()
   [2, 3]
   sage: P2 = Kgood.prime_above(2)
   sage: N2 = Egamma.conductor_exponent(P2); N2
   12 if ('a', 'b') == (1, 0) mod 2
   14 if ('a', 'b') == (1, 1) mod 2
   8  if ('a', 'b') == (0, 3), (2, 3) mod 4
   0  if ('a', 'b') is 1 of 4 possibilities mod 8
   4  if ('a', 'b') is 1 of 4 possibilities mod 8
   sage: N2[3]
   (0, The condition that ('a', 'b') == (0, 1), (2, 5), (4, 1), (6, 5) mod 8)
   sage: P3 = Kgood.prime_above(3)
   sage: N3 = Egamma.conductor_exponent(P3); N3
   0 if ('a', 'b') is 1 of 6 possibilities mod 3
   2 if ('a', 'b') == (0, 1), (0, 2) mod 3
   sage: N3[0]
   (0,
    The condition that ('a', 'b') == (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2) mod 3)

To summarize we now have:
 - additive reduction at the prime above 2, unless :math:`a \equiv 0`
   (mod 2) and :math:`b \equiv a^2 + 1` (mod 8) in which case the
   reduction is good.
 - additive reduction at the prime above 3, unless :math:`3 \nmid a`
   in which case the reduction is good.
 - at any other prime :math:`P`, good reduction if :math:`P \nmid
   \Delta` and multiplicative reduction otherwise.

Now we determine minimal models at the primes above 2 and 3, to
determine minimal models everywhere.

::

   sage: Egamma2 = Egamma.minimal_model(P2); Egamma2
   Elliptic Curve defined by y^2 = x^3 + ((6*lu0+12)*a)*x^2 + ((81/2*lu0+162)*a^2+(9*lu0+27)*b)*x + ((81*lu0+270)*a^3+(135*lu0+486)*a*b) over Multivariate Polynomial Ring in a, b over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu                                                                                                                                                                                                                 if ('a', 'b') is 1 of 44 possibilities mod 8
   Elliptic Curve defined by y^2 + ((-1/2*lu0-1))*x*y + ((24*lu0+84)*a+(5/2*lu0+8))*y = x^3 + ((-18*lu0-60)*a+(-lu0-4))*x^2 + ((5535/8*lu0+4797/2)*a^2+(99*lu0+342)*a+(117/4*lu0+405/4)*b+(19/4*lu0+67/4))*x + ((-70227/8*lu0-121635/4)*a^3+(-4599/2*lu0-63729/8)*a^2+(-7155/8*lu0-12393/4)*a*b+(-525/2*lu0-909)*a+(-873/16*lu0-189)*b+(-175/16*lu0-38)) over Multivariate Polynomial Ring in a, b over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu if ('a', 'b') is 1 of 4 possibilities mod 8
   sage: Egamma2[1]
   (Elliptic Curve defined by y^2 + ((-1/2*lu0-1))*x*y + ((24*lu0+84)*a+(5/2*lu0+8))*y = x^3 + ((-18*lu0-60)*a+(-lu0-4))*x^2 + ((5535/8*lu0+4797/2)*a^2+(99*lu0+342)*a+(117/4*lu0+405/4)*b+(19/4*lu0+67/4))*x + ((-70227/8*lu0-121635/4)*a^3+(-4599/2*lu0-63729/8)*a^2+(-7155/8*lu0-12393/4)*a*b+(-525/2*lu0-909)*a+(-873/16*lu0-189)*b+(-175/16*lu0-38)) over Multivariate Polynomial Ring in a, b over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu,
    The condition that ('a', 'b') == (0, 1), (2, 5), (4, 1), (6, 5) mod 8)
   sage: from modular_method.diophantine_equations.conditions import apply_to_conditional_value
   sage: apply_to_conditional_value(
   ....:     lambda Egamma_: (Egamma.discriminant() / Egamma_.discriminant()).factor(),
   ....:     Egamma2)
   (-1053780*lu0 + 3650401) * (1/2*lu0 + 1)^12              if ('a', 'b') is 1 of 44 possibilities mod 8
   (-7693439131560*lu0 + 26650854921601) * (1/2*lu0 + 1)^24 if ('a', 'b') is 1 of 4 possibilities mod 8
   sage: Egamma3 = Egamma.minimal_model(P3); Egamma3
   Elliptic Curve defined by y^2 = x^3 + ((-1/2*lu0+1)*a)*x^2 + ((-1/8*lu0+1/2)*a^2+(1/4*lu0-3/4)*b)*x over Multivariate Polynomial Ring in a, b over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu
   sage: (Egamma.discriminant() / Egamma3.discriminant()).factor()
   (-1053780*lu0 + 3650401) * 1/2*lu0^12 * (1/2*lu0 + 1)^24
   sage: P2, P3
   (Fractional ideal (1/2*lu0 + 1), Fractional ideal (1/2*lu0))

We thus see that to make a global minimal model of `Egamma` we need to
scale the curve by a generator of `P2` times a generator of `P3`, and
once again with a generator of `P2` when :math:`a \equiv 0` (mod 2)
and :math:`b \equiv a^2 + 1` (mod 8). This allows us to define a
global minimal model `Egood` for the curve.

::
   
   sage: pi2 = P2.gens_reduced()[0]
   sage: pi3 = P3.gens_reduced()[0]
   sage: Egood_ = FreyCurve(Egamma.scale_curve((pi2*pi3)^(-1)), condition=Egamma._condition)
   sage: Egood = Egood_.minimal_model(P2)

We turn these minimal models into Q-curves by constructing the
corresponding isogenies from the isogenies of `Egamma` combined with
the isomorphisms from `Egood` to `Egamma`

::

   sage: Ggood.<sgood> = Kgood.galois_group()
   sage: def make_Qcurve(E, con):
   ....:     u = sqrt(Kgood(E.c4().parent()((E.c4() * Egamma.c6()) / (Egamma.c4() * E.c6()))))
   ....:     r = (u^2*E.b2() - Egamma.b2()) / 12
   ....:     su = sgood(u)
   ....:     sr = r.change_ring(sgood.as_hom())
   ....:     F = Egamma.isogeny_x_map(sgood)
   ....:     l = Kgood(Egamma.isogeny_scalar(sgood))
   ....:     x = F.parent().gen()
   ....:     Fnew = (F(u^2*x + r) - sr) / su^2
   ....:     lnew = l * u^(-1) * su
   ....:     return FreyQcurve(E, isogenies={sgood: (Fnew, lnew)}, condition=Egamma._condition & con)
   sage: Egood = apply_to_conditional_value(make_Qcurve, Egood, use_condition=True)
   sage: apply_to_conditional_value(
   ....:     lambda E_: E_.minimal_model(P2).a_invariants() == E_.a_invariants(),
   ....:     Egood)
   True
   sage: apply_to_conditional_value(
   ....:     lambda E_: E_.minimal_model(P3).a_invariants() == E_.a_invariants(),
   ....:     Egood)
   True
   sage: apply_to_conditional_value(
   ....:     lambda E_: set(P.smallest_integer()
   ....:         for P, _ in Kgood.ideal(E_.c4().resultant(E_.discriminant())(1, 1)).factor()),
   ....:     Egood)
   {2, 3}

We thus obtained global minimal models for each possible case. Next we
check what the x-maps associated with each isogeny are, to see for
which primes `P` the reduction of the isogeny modulo `P` is separable.

::

   sage: apply_to_conditional_value(lambda E_: E_.isogeny_x_map(Ggood(1)), Egood)
   x
   sage: a, b = Egamma.parameters()
   sage: F0 = Egood[0][0].isogeny_x_map(sgood)
   sage: x = F0.parent().gen()
   sage: F0 == ((26 + 15*sqrt3)*x^2 - (10 + 6*sqrt3)*a*x + (a^2 + b*sqrt3)) / (2*x)
   True
   sage: F1 = Egood[1][0].isogeny_x_map(sgood)
   sage: x = F1.parent().gen()
   sage: F1 == (((2 + sqrt3)*x^2 +
   ....:         ((3 + 3*sqrt3)*a - (5 + 5*sqrt3))*x +
   ....:         17*(a/2)^2 - (25 + 3*sqrt3)*(a/2) + sqrt3*(b - 1)/4 + (6 + 4*sqrt3)) /
   ....:        (2*x + (-2 + 2*sqrt3)*a - (2 + sqrt3)))
   True

It is easy to see that for all primes `P` of characteristic
:math:`\neq 2` these x-maps would reduce to maps with non-zero
derivative modulo `P`, hence the reduction of the corresponding
isogenies is separable by Proposition 2.10.2. The same is not true
when `P` is the prime above 2. In that case `F0` does not reduce at
all, and `F1` reduces to the inseparable map `x^2`.

Since we can not solve the case :math:`p = 2` we will from now on work
with the model `Egood[0][0]`.

::

   sage: Egood = Egood[0][0]

We now apply Theorem 2.10.1 and Proposition 2.10.3 to compute the
trace of Frobenius for each odd prime number :math:`p` and the default
splitting map. To ease the process we define some intermediary
functions. First a function that computes the condition for which we
have good reduction at an odd prime :math:`p`. 

::

   sage: from modular_method.diophantine_equations.conditions import ConditionalValue
   sage: def good_con(p):
   ....:     has_good_red = Egood.has_good_reduction(Kgood.prime_above(p))
   ....:     if isinstance(has_good_red, ConditionalValue):
   ....:         return next(con for val, con in has_good_red if val)
   ....:     elif has_good_red:
   ....:         return Egood._condition
   ....:     else:
   ....:         a, b = Egood.parameters()
   ....:         a, b = a.change_ring(QQ), b.change_ring(QQ)
   ....:         return ~CoprimeCondition([a, b], 0)

Next a function that computes the possible reductions
:math:`\tilde{E}` together with the reduction of the x-map `F` at a
prime above :math:`p`. For the reduction of the map `F` we use the
formula for the numerator and denominator of `F0` found before.

::

   sage: from modular_method.padics.pAdic_base import pAdicBase
   sage: from modular_method.padics.pAdic_tree import pAdicNode, pAdicTree
   sage: from modular_method.diophantine_equations.conditions import TreeCondition
   sage: def reduction_data(p, a, b):
   ....:      P = Kgood.prime_above(p)
   ....:      Ered = Egood.specialize((a, b)).reduction(P)
   ....:      _.<x> = Kgood[]
   ....:      FP = P.residue_field()
   ....:      Fnum = ((26 + 15*sqrt3)*x^2 - (10 + 6*sqrt3)*a*x + (a^2 + b*sqrt3)).change_ring(FP)
   ....:      Fden = (2*x).change_ring(FP)
   ....:      return Ered, Fnum / Fden
   sage: def possible_reductions(p):
   ....:      P = Kgood.prime_above(p)
   ....:      pAdics = pAdicBase(QQ, p)
   ....:      T = good_con(p).pAdic_tree(pAdics=pAdics)
   ....:      node_ls, root = T.nodes_at_level(1)
   ....:      result = ConditionalValue([(
   ....:          reduction_data(p, *node.representative()),
   ....:          TreeCondition(pAdicTree([a, b], root=node.sub_tree())))
   ....:          for node in node_ls])
   ....:      return apply_to_conditional_value(lambda E: E, result)

We turn the computation of the trace into functions as well, so they
work well with `apply_to_conditional_value`. First the case when
`\sigma` acts trivial on `Kgood`. Here the trace directly follows as
in Example 2.10.4.

::

   sage: def trace1(p, Ered, Fred):
   ....:     return 1 + p - Ered.count_points()

Next the case when `\sigma` acts non-trivially on `Kgood`, in which
case we need to compute the quantities in Proposition 2.10.3. The
cardinality of the set computed in Proposition 2.10.3 is denoted below
by `m`

::

   sage: def trace2(p, Ered, Fred):
   ....:     beta = Egood.splitting_map()
   ....:     x = Fred.numerator().parent().gen()
   ....:     FP = Ered.base()
   ....:     f1 = (Fred - x^p).numerator()
   ....:     R = 4*x^3 + Ered.b2()*x^2 + 2*Ered.b4()*x + Ered.b6()
   ....:     l = FP(Egood.isogeny_scalar(sgood))
   ....:     f2 = (l*R^((p+1)/2) - Fred.derivative(x)*R).numerator()
   ....:     m = 1 + 2*gcd(f1, f2).radical().degree() - gcd(f1, R).radical().degree()
   ....:     aEp = 2 + p - m
   ....:     return beta(sgood)^(-1) * aEp
   
With these functions we can now do a similar loop as for the non-Frey
curve case, where some parts are replaced by the functions above.

::

   sage: Lbeta.<sqrtm2> = QuadraticField(-2)
   sage: for p in prime_range(3, 30):
   ....:     red_data = possible_reductions(p)
   ....:     P = Kgood.prime_above(p)
   ....:     FP = P.residue_field()
   ....:     if FP.degree() == 1:
   ....:         # The case sigma in G_K
   ....:         trace = apply_to_conditional_value(lambda data: trace1(p, *data), red_data)
   ....:         print(p, "ramifies/splits, trace:")
   ....:         print(trace)
   ....:     if len(Kgood.primes_above(p)) == 1:
   ....:         # The case sigma not in G_K
   ....:         trace = apply_to_conditional_value(lambda data: trace2(p, *data), red_data)
   ....:         print(p, "ramifies/inert, trace:")
   ....:         print(trace)
   3 ramifies/splits, trace:
   -2 if ('a', 'b') is 1 of 3 possibilities mod 3
   2  if ('a', 'b') is 1 of 3 possibilities mod 3
   3 ramifies/inert, trace:
   -2*zeta4a0 if ('a', 'b') is 1 of 3 possibilities mod 3
   2*zeta4a0  if ('a', 'b') is 1 of 3 possibilities mod 3
   5 ramifies/inert, trace:
   0          if ('a', 'b') is 1 of 4 possibilities mod 5
   zeta4a0    if ('a', 'b') is 1 of 5 possibilities mod 5
   -3*zeta4a0 if ('a', 'b') == (0, 2) mod 5
   2*zeta4a0  if ('a', 'b') is 1 of 4 possibilities mod 5
   3*zeta4a0  if ('a', 'b') == (0, 3) mod 5
   -2*zeta4a0 if ('a', 'b') is 1 of 4 possibilities mod 5
   -zeta4a0   if ('a', 'b') is 1 of 5 possibilities mod 5
   7 ramifies/inert, trace:
   0          if ('a', 'b') is 1 of 12 possibilities mod 7
   -3*zeta4a0 if ('a', 'b') is 1 of 6 possibilities mod 7
   zeta4a0    if ('a', 'b') is 1 of 6 possibilities mod 7
   -2*zeta4a0 if ('a', 'b') is 1 of 6 possibilities mod 7
   2*zeta4a0  if ('a', 'b') is 1 of 6 possibilities mod 7
   -zeta4a0   if ('a', 'b') is 1 of 6 possibilities mod 7
   3*zeta4a0  if ('a', 'b') is 1 of 6 possibilities mod 7
   11 ramifies/splits, trace:
   6  if ('a', 'b') is 1 of 5 possibilities mod 11
   -6 if ('a', 'b') is 1 of 5 possibilities mod 11
   0  if ('a', 'b') is 1 of 30 possibilities mod 11
   -4 if ('a', 'b') is 1 of 20 possibilities mod 11
   -2 if ('a', 'b') is 1 of 10 possibilities mod 11
   2  if ('a', 'b') is 1 of 10 possibilities mod 11
   4  if ('a', 'b') is 1 of 20 possibilities mod 11
   13 ramifies/splits, trace:
   0  if ('a', 'b') is 1 of 12 possibilities mod 13
   -4 if ('a', 'b') is 1 of 15 possibilities mod 13
   2  if ('a', 'b') is 1 of 36 possibilities mod 13
   -2 if ('a', 'b') is 1 of 36 possibilities mod 13
   -6 if ('a', 'b') is 1 of 15 possibilities mod 13
   6  if ('a', 'b') is 1 of 15 possibilities mod 13
   4  if ('a', 'b') is 1 of 15 possibilities mod 13
   17 ramifies/inert, trace:
   4*zeta4a0  if ('a', 'b') is 1 of 24 possibilities mod 17
   -4*zeta4a0 if ('a', 'b') is 1 of 24 possibilities mod 17
   -5*zeta4a0 if ('a', 'b') is 1 of 20 possibilities mod 17
   -2*zeta4a0 if ('a', 'b') is 1 of 32 possibilities mod 17
   0          if ('a', 'b') is 1 of 32 possibilities mod 17
   5*zeta4a0  if ('a', 'b') is 1 of 20 possibilities mod 17
   2*zeta4a0  if ('a', 'b') is 1 of 32 possibilities mod 17
   3*zeta4a0  if ('a', 'b') is 1 of 20 possibilities mod 17
   -zeta4a0   if ('a', 'b') is 1 of 32 possibilities mod 17
   zeta4a0    if ('a', 'b') is 1 of 32 possibilities mod 17
   -3*zeta4a0 if ('a', 'b') is 1 of 20 possibilities mod 17
   19 ramifies/inert, trace:
   -6*zeta4a0 if ('a', 'b') is 1 of 9 possibilities mod 19
   6*zeta4a0  if ('a', 'b') is 1 of 9 possibilities mod 19
   0          if ('a', 'b') is 1 of 54 possibilities mod 19
   3*zeta4a0  if ('a', 'b') is 1 of 54 possibilities mod 19
   -zeta4a0   if ('a', 'b') is 1 of 18 possibilities mod 19
   2*zeta4a0  if ('a', 'b') is 1 of 36 possibilities mod 19
   -4*zeta4a0 if ('a', 'b') is 1 of 18 possibilities mod 19
   -5*zeta4a0 if ('a', 'b') is 1 of 18 possibilities mod 19
   -3*zeta4a0 if ('a', 'b') is 1 of 54 possibilities mod 19
   5*zeta4a0  if ('a', 'b') is 1 of 18 possibilities mod 19
   4*zeta4a0  if ('a', 'b') is 1 of 18 possibilities mod 19
   -2*zeta4a0 if ('a', 'b') is 1 of 36 possibilities mod 19
   zeta4a0    if ('a', 'b') is 1 of 18 possibilities mod 19
   23 ramifies/splits, trace:
   0  if ('a', 'b') is 1 of 132 possibilities mod 23
   -8 if ('a', 'b') is 1 of 44 possibilities mod 23
   -2 if ('a', 'b') is 1 of 22 possibilities mod 23
   8  if ('a', 'b') is 1 of 44 possibilities mod 23
   -4 if ('a', 'b') is 1 of 66 possibilities mod 23
   -6 if ('a', 'b') is 1 of 44 possibilities mod 23
   4  if ('a', 'b') is 1 of 66 possibilities mod 23
   6  if ('a', 'b') is 1 of 44 possibilities mod 23
   2  if ('a', 'b') is 1 of 22 possibilities mod 23
   29 ramifies/inert, trace:
   0          if ('a', 'b') is 1 of 28 possibilities mod 29
   7*zeta4a0  if ('a', 'b') is 1 of 35 possibilities mod 29
   5*zeta4a0  if ('a', 'b') is 1 of 56 possibilities mod 29
   -6*zeta4a0 if ('a', 'b') is 1 of 28 possibilities mod 29
   -4*zeta4a0 if ('a', 'b') is 1 of 56 possibilities mod 29
   -5*zeta4a0 if ('a', 'b') is 1 of 56 possibilities mod 29
   zeta4a0    if ('a', 'b') is 1 of 56 possibilities mod 29
   2*zeta4a0  if ('a', 'b') is 1 of 112 possibilities mod 29
   -2*zeta4a0 if ('a', 'b') is 1 of 112 possibilities mod 29
   -3*zeta4a0 if ('a', 'b') is 1 of 63 possibilities mod 29
   3*zeta4a0  if ('a', 'b') is 1 of 63 possibilities mod 29
   4*zeta4a0  if ('a', 'b') is 1 of 56 possibilities mod 29
   -7*zeta4a0 if ('a', 'b') is 1 of 35 possibilities mod 29
   -zeta4a0   if ('a', 'b') is 1 of 56 possibilities mod 29
   6*zeta4a0  if ('a', 'b') is 1 of 28 possibilities mod 29

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
