========================================
 Sur l'équation :math:`a^3 + b^3 = c^p`
========================================

We run some computations for results in the article "Sur l'équation
:math:`a^3 + b^3 = c^p` written by Alain Kraus and published in
Experimental Mathematics, Volume 7 (1998), Issue 1.

.. linkall

::

   sage: from modular_method import *

The article considers the equation :math:`a^3 + b^3 = c^p` with
:math:`a, b, c` pairwise coprime integers, :math:`p \ge 3` a prime
number, and :math:`a b c \ne 0`. We will use :math:`a` and :math:`b`
as parameters in the code and write `cp` for :math:`c^p`.

::

   sage: R.<a, b> = QQ[]
   sage: cp = a^3 + b^3

In section 4 the article introduces a Frey curve to solve this
Diophantine equation. Before they do they also make the additional
assumptions that :math:`a c` is even and :math:`p \ge 5`. We introduce
this Frey curve here with the corresponding condition.

::

   sage: condition = (CoprimeCondition([a, b]) &
   ....:              CongruenceCondition(a*cp, 2) &
   ....:              PowerCondition(cp, 5))
   sage: Eab = FreyCurve([0, 0, 0, 3*a*b, b^3 - a^3], condition=condition)

We check that the :math:`c_4`, :math:`c_6` and :math:`\Delta`
invariants agree with the article.

::

   sage: Eab.c4() == - 2^4 * 3^2 * a * b
   True
   sage: Eab.c6() == 2^5 * 3^3 * (a^3 - b^3)
   True
   sage: Eab.discriminant() == -2^4 * 3^3 * cp^2
   True

Next we compute the conductor of this Frey curve.

::

   sage: N = Eab.conductor(); N
   Warning: Assuming that a and b are coprime.
   2^n0*9*Rad_P( (-432) * (a + b)^2 * (a^2 - a*b + b^2)^2 )
    where 
   n0 = 2 if ('a', 'b') == (0, 1) mod 4
        4 if ('a', 'b') is 1 of 136 possibilities mod 32
        3 if ('a', 'b') == (2, 1) mod 4
        1 if ('a', 'b') is 1 of 8 possibilities mod 32
   sage: cp.factor()
   (a + b) * (a^2 - a*b + b^2)
   sage: N.left().left().right()[3][1]
   The condition that ('a', 'b') == (1, 31), (5, 27), (9, 23), (13, 19), (17, 15), (21, 11), (25, 7), (29, 3) mod 32

It can easily be checked that this agrees with the result given in
Lemme 4.1 of the article. The article discusses in Section 5 why the
left hand side will be the level of associated newforms of this
elliptic curve. We will compute the newforms directly, taking into
account the additional assumptions in Section~6.

::

   sage: extra = ((CongruenceCondition(cp, 2) & CongruenceCondition(b + 1, 4)) |
   ....:          (~CongruenceCondition(cp, 2) & CongruenceCondition(b - 1, 4)))
   sage: nfs = Eab.newform_candidates(condition=condition & extra); nfs   
   Warning: The bad primes chosen by default only take into account primes of additive reduction.
   [q + O(q^6)]         if ('a', 'b') == (0, 1) mod 4
   [q + 2*q^5 + O(q^6)] if ('a', 'b') == (2, 1) mod 4
   []                   if ('a', 'b') is 1 of 8 possibilities mod 32

We compute the elliptic curves corresponding to these rational newforms.

::

   sage: nfs[0][0][0]._f.abelian_variety().elliptic_curve()
   Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
   sage: nfs[1][0][0]._f.abelian_variety().elliptic_curve()
   Elliptic Curve defined by y^2 = x^3 + 6*x - 7 over Rational Field

Note that these correspond to the curves :math:`E'` introduced in
Section 6 and :math:`E` introduced in Section 2B. Therefore these
computations corroborate Proposition 6.3.

What is not shown explicitly in the article, but what we will show
here is that normal elimination by traces seems to yield no results.

::

   sage: nfs = eliminate_by_traces(Eab, nfs, primes=prime_range(5, 50),
   ....:                           condition=CoprimeCondition([a, b])); nfs
   [(q + O(q^6), 0)]         if ('a', 'b') == (0, 1) mod 4
   [(q + 2*q^5 + O(q^6), 0)] if ('a', 'b') == (2, 1) mod 4
   []                        if ('a', 'b') is 1 of 8 possibilities mod 32

In Thèorème 6.1 in the article it is shown that we are in the case
where :math:`a \equiv 2` modulo 4, :math:`2 \nmid c` and :math:`3 \mid
c`. This implies that :math:`3 (a + b) = c_1^p`, :math:`a^2 - a b +
b^2 = 3 c_2^p` and :math:`c = c_1 c_2` with :math:`3 \mid c_1` and
:math:`3 \nmid c_2`. We take this further by factoring over the
cyclotomic field :math:`\Q(\zeta_3)` showing that :math:`a + \zeta_3 b
= (1 + 2 \zeta_3) \gamma^p` with :math:`c_2 = \gamma \bar{\gamma}`.

The article uses this information to eliminate every prime between 16
and 10000. The lower bound on the prime :math:`p` seems to arise from
the irreducibility result used in Lemme 5.2. We will assume that the
representation associated with `Eab` is irreducible for all primes
:math:`p \ge 5` and do the kraus method for all primes :math:`p`
between 4 and 50, as this already takes a while for the framework. The
article can go beyond this by going from two to one parameter. The
framework can not do this automatically so we stick to the part it can
do automatically.

::

   sage: nf = nfs[1][0][0][0]; nf
   q + 2*q^5 + O(q^6)
   sage: K.<zeta3> = CyclotomicField(3)
   sage: c1p = 3*(a + b)
   sage: c2p = (a + zeta3*b) / (1 + 2*zeta3)
   sage: all(len(kraus_method(Eab, [(nf, p)], p, [c1p, c2p], primes=30*p+2,
   ....:                      condition=CoprimeCondition([a, b]))) == 0
   ....:     for p in prime_range(5, 50))
   True
