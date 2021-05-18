=================================================================================
 A Multi-Frey Approach to Some Multi-Parameter families of Diophantine Equations
=================================================================================

We run the computations for some results in the article "A Multi-Frey
Approach to Some Multi-Parameter families of Diophantine Equations"
written by Yann Bugeaud, Maurice Mignotte, and Samir Siksek, and
published in Canadian Journal of Mathematics, volume 60 (2008),
issue 3. The article can be found at `Samir Siksek's website`_

.. _Samir Siksek's website: https://homepages.warwick.ac.uk/staff/S.Siksek/papers/tripexp13.pdf
.. linkall

The article considers Diophantine equations of the form :math:`\alpha
x^p - 2^r \beta y^p = 1` with :math:`\alpha, \beta, x` odd and
:math:`y` non-zero. They furthermore assume that :math:`alpha` and
:math:`beta` are :math:`p` th power free, and that :math:`p \ge
7`. Throughout the paper they use :math:`\psi = 2^r \beta y^p` as a
parameter for Frey curves which they assume to be unequal to
:math:`\pm 2`.

::

   sage: from modular_method import *
   sage: R.<psi> = QQ[]
   sage: con_base = CongruenceCondition(psi, 2)

Next they introduce a variety of Frey curves with parameter :math:`\psi`

::

   sage: F1 = FreyCurve([0, 2*psi + 1, 0, psi^2 + psi, 0], condition=con_base)
   sage: F2 = FreyCurve([0, -(2*psi + 1), 0, psi^2 + psi, 0], condition=con_base)
   sage: G1 = FreyCurve([1, 0, 0, -psi/64, 0], condition=CongruenceCondition(psi, 64))
   sage: G2 = FreyCurve([0, 1, 0, -psi/4, 0], condition=CongruenceCondition(psi, 4))
   sage: G3 = FreyCurve([0, -1, 0, -psi/4, 0], condition=CongruenceCondition(psi, 4))
   sage: G4 = FreyCurve([0, 2, 0, -psi, 0], condition=con_base)
   sage: H1 = FreyCurve([3, 0, -psi, 0, 0], condition=con_base)
   sage: H2 = FreyCurve([3, 0, psi + 1, 0, 0], condition=con_base)
   sage: H3 = FreyCurve([-3, 0, psi, 0, 0], condition=con_base)

In the tables 1 through 4 these Frey curves are assigned to specific
cases for the variables. Note that the cases (iii) and (iv) in table 4
do not appear as they imply :math:`3 \mid \alpha`. We will define the
elliptic curves through conditional values.

::

   sage: from modular_method.diophantine_equations.conditions import ConditionalValue
   sage: F = ConditionalValue([(F1, CongruenceCondition(psi, 4)),
   ....:                       (F2, ~CongruenceCondition(psi, 4) &
   ....:                            CongruenceCondition(psi, 2))])
   sage: G = ConditionalValue([(G1, CongruenceCondition(psi, 64)),
   ....:                       (G2, ~CongruenceCondition(psi, 64) &
   ....:                            (CongruenceCondition(psi, 8) |
   ....:                             CongruenceCondition(psi - 4, 16))),
   ....:                       (G3, CongruenceCondition(psi + 4, 16)),
   ....:                       (G4, ~CongruenceCondition(psi, 4) &
   ....:                            CongruenceCondition(psi, 2))])
   sage: H = ConditionalValue([(H1, CongruenceCondition(psi, 3) |
   ....:                            CongruenceCondition(psi - 4, 9)),
   ....:                       (H2, CongruenceCondition(psi + 1, 3)),
   ....:                       (H3, CongruenceCondition(psi - 1, 9) |
   ....:                            CongruenceCondition(psi - 7, 9))])

We now confirm the result about the levels of corresponding newforms
claimed in Proposition 3.2. First we determine the primes at which
these curves can have additive reduction.

::

   sage: from modular_method.diophantine_equations.conditions import apply_to_conditional_value
   sage: apply_to_conditional_value(lambda Fi : Fi.primes_of_possible_additive_reduction(), F)
   [2]
   sage: apply_to_conditional_value(lambda Gi : Gi.primes_of_possible_additive_reduction(), G)   
   []  if psi == 0 mod 64
   [2] if psi ~= 0 mod 64 and psi == 0 mod 8 or psi - 4 == 0 mod 16 or psi + 4 == 0 mod 16 or psi ~= 0 mod 4 and psi == 0 mod 2
   sage: apply_to_conditional_value(lambda Hi : Hi.primes_of_possible_additive_reduction(), H)
   [3]

This shows that the conductors are a power of 2, or a power of 3 for
:math:`H`, times the radical of the discriminant outside these
primes. We print the discriminants here, to show that they are indeed
all of the form :math:`\psi^m (\psi + 1)^n` times some power of 2 or 3
as claimed in the article.

::
   
   sage: apply_to_conditional_value(lambda Fi : Fi.discriminant().factor(), F)
   (16) * psi^2 * (psi + 1)^2
   sage: apply_to_conditional_value(lambda Gi : Gi.discriminant().factor(), G)   
   (1/4096) * (psi + 1) * psi^2 if psi == 0 mod 64
   (psi + 1) * psi^2            if psi ~= 0 mod 64 and psi == 0 mod 8 or psi - 4 == 0 mod 16 or psi + 4 == 0 mod 16
   (64) * (psi + 1) * psi^2     if psi ~= 0 mod 4 and psi == 0 mod 2
   sage: apply_to_conditional_value(lambda Hi : Hi.discriminant().factor(), H)
   (-27) * (psi + 1) * psi^3 if psi == 0 mod 3 or psi - 4 == 0 mod 9 or psi - 1 == 0 mod 9 or psi - 7 == 0 mod 9
   (-27) * psi * (psi + 1)^3 if psi + 1 == 0 mod 3

Next we compute the conductor exponent at the corresponding possible
additive prime: 2 for :math:`F, G` and 3 for :math:`H`.
   
::

   sage: e1 = apply_to_conditional_value((lambda Fi, con : Fi.conductor_exponent(2, condition=con)),
   ....:                                 F, use_condition=True, default_condition=con_base)
   sage: L1 = 2^e1; L1
   2^n0
    where 
   n0 = 3 if psi == 4, 8, 12 mod 16
        0 if psi == 16 mod 32
        1 if psi == 0 mod 32
        5 if psi ~= 0 mod 4 and psi == 0 mod 2
   sage: e2 = apply_to_conditional_value((lambda Gi, con : Gi.conductor_exponent(2, condition=con)),
   ....:                                 G, use_condition=True, default_condition=con_base)
   sage: L2 = 2^e2; L2
   2^n0
    where 
   n0 = 0 if psi == 64 mod 128
        1 if psi == 0 mod 128
        2 if psi == 4 mod 16
        5 if psi == 8 mod 16
        3 if psi == 16, 32, 48 mod 64 or psi + 4 == 0 mod 16
        7 if psi ~= 0 mod 4 and psi == 0 mod 2
   sage: e3 = apply_to_conditional_value((lambda Hi, con : Hi.conductor_exponent(3, condition=con)),
   ....:                                 H, use_condition=True, default_condition=con_base)
   sage: L3 = 2^e3; L3
   2^n0
    where 
   n0 = 2 if psi == 4 mod 9
        4 if psi == 2, 3, 5, 6 mod 9
        3 if psi == 8, 9, 17, 18 mod 27 or psi - 1 == 0 mod 9 or psi - 7 == 0 mod 9
        0 if psi == 26, 27, 53, 54 mod 81
        1 if psi == 0, 80 mod 81

These agree with the values of :math:`L_1, L_2, L_3` in the tables 1
through 4, except for case (II) in table 1, and case (i) in tables 2
and 3 when :math:`3 \mid y` or :math:`3 \mid x`. These seem to be
mistakes in the original article.

Now we do the newform computations for the curves obtained in the
first example of Section 6. First we determine all the appropriate
cases

::

   sage: from modular_method.diophantine_equations.conditions import conditional_product
   sage: N1 = (5*L1).value()
   sage: N2 = (5*L2).value()
   sage: cases = conditional_product(F, G, N1, N2)

We write each case a little bit simpler by replacing with tree
conditions and skipping the empty ones.

::

   sage: from modular_method.padics.pAdic_base import pAdicBase
   sage: pAdics = pAdicBase(QQ, 2)
   sage: cases = [(val, TreeCondition(con.pAdic_tree(pAdics=pAdics))) for val, con in cases
   ....:          if not con.pAdic_tree(pAdics=pAdics).is_empty()]

Now we do all the cases at once

::

   sage: for (F_, G_, N1, N2), con in cases:
   ....:     nfs = [(f, g) for f in get_newforms(N1) for g in get_newforms(N2)]
   ....:     nfs = eliminate_by_traces((F_, G_), nfs, primes=[3, 7, 11, 13])
   ....:     print(con)
   ....:     for f, g, n in nfs:
   ....:         print("    ", f, "|", g, "|", n)
   The condition that psi == 64 mod 128
   The condition that psi == 0 mod 128
   The condition that psi == 4 mod 16
        q + q^5 + O(q^6) | q - 2*q^3 - q^5 + O(q^6) | 0
   The condition that psi == 8 mod 16
        q + q^5 + O(q^6) | q - 2*q^3 - q^5 + O(q^6) | 56
        q + q^5 + O(q^6) | q + 2*q^3 - q^5 + O(q^6) | 12
        q + q^5 + O(q^6) | q + 1/2*a2*q^3 + q^5 + O(q^6) | 4
   The condition that psi == 16 mod 32
   The condition that psi == 32 mod 64
   The condition that psi == 12 mod 16
        q + q^5 + O(q^6) | q + q^5 + O(q^6) | 4
   The condition that psi == 2 mod 4
        q - 2*q^3 - q^5 + O(q^6) | q - 2*q^3 - q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q - 2*q^3 + q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q - q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q - q^5 + O(q^6) | 6
        q - 2*q^3 - q^5 + O(q^6) | q + q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q + q^5 + O(q^6) | 6
        q - 2*q^3 - q^5 + O(q^6) | q + 2*q^3 - q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q + 2*q^3 + q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q - 1/2*a8*q^3 - q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q + 1/2*a9*q^3 + q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q + 1/2*a10*q^3 - q^5 + O(q^6) | 2
        q - 2*q^3 - q^5 + O(q^6) | q + 1/2*a11*q^3 + q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q - 2*q^3 - q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q - 2*q^3 + q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q - q^5 + O(q^6) | 6
        q + 2*q^3 - q^5 + O(q^6) | q - q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q + q^5 + O(q^6) | 6
        q + 2*q^3 - q^5 + O(q^6) | q + q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q + 2*q^3 - q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q + 2*q^3 + q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q - 1/2*a8*q^3 - q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q + 1/2*a9*q^3 + q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q + 1/2*a10*q^3 - q^5 + O(q^6) | 2
        q + 2*q^3 - q^5 + O(q^6) | q + 1/2*a11*q^3 + q^5 + O(q^6) | 2
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q - 2*q^3 - q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q - 2*q^3 + q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q - q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q - q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q + q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q + q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q + 2*q^3 - q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q + 2*q^3 + q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q - 1/2*a8*q^3 - q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q + 1/2*a9*q^3 + q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q + 1/2*a10*q^3 - q^5 + O(q^6) | 4
        q + 1/2*a2*q^3 + q^5 + O(q^6) | q + 1/2*a11*q^3 + q^5 + O(q^6) | 4

We see that indeed only for :math:`\psi \equiv 4` modulo 16 (case VII
in the article) we can not eliminate all newforms for primes :math:`p
\ge 7`.

We now do the same to the second example of Section 6, noting that in
this case we should include :math:`H`. Note that the level
corresponding to :math:`H` partially depends on whether :math:`r >
0`. Since we do not have :math:`r` as a variable we add a
TextCondition for this.

::

   sage: from modular_method.diophantine_equations.conditions import TextCondition
   sage: rgt0 = ConditionalValue([(1, TextCondition('r > 0')),
   ....:                          (0, TextCondition('r = 0'))])
   sage: N1 = (3*L1).value()
   sage: N2 = (3*L2).value()
   sage: N3 = (L3 * 2^rgt0).value()
   sage: cases = conditional_product(F, G, H, N1, N2, N3)

Again we will write each case a bit simpler using tree
conditions. Note that to reintroduce the TextCondition we look at
whether N3 is divisible by 2.

::

   sage: pAdics2 = pAdicBase(QQ, 2)
   sage: pAdics3 = pAdicBase(QQ, 3)
   sage: cases = [((F_, G_, H_, N1, N2, N3),
   ....:           (TreeCondition(con.pAdic_tree(pAdics=pAdics2)) &
   ....:            TreeCondition(con.pAdic_tree(pAdics=pAdics3)) &
   ....:            (TextCondition('r > 0') if 2.divides(N3)
   ....:             else TextCondition('r = 0'))))
   ....:          for (F_, G_, H_, N1, N2, N3), con in cases
   ....:          if not con.pAdic_tree(pAdics=pAdics2).is_empty()
   ....:          and not con.pAdic_tree(pAdics=pAdics3).is_empty()]

Again we now do all the cases at once and obtain similar results as in
the article

::

   sage: for (F_, G_, H_, N1, N2, N3), con in cases:
   ....:     nfs = [(f, g, h) for f in get_newforms(N1)
   ....:            for g in get_newforms(N2) for h in get_newforms(N3)]
   ....:     nfs = eliminate_by_traces((F_, G_, H_), nfs, primes=[5, 7, 11])
   ....:     if len(nfs) > 0:
   ....:         print(con)
   ....:         for f, g, h, n in nfs:
   ....:              print("    ", f, "|", g, "|", h, "|", n)
   The condition that psi == 8 mod 16 and the condition that psi == 3, 6 mod 9 and r > 0
        q - q^3 - 2*q^5 + O(q^6) | q - q^3 + 2*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 4
        q - q^3 - 2*q^5 + O(q^6) | q + q^3 + 2*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 4
   The condition that psi == 8 mod 16 and the condition that psi == 2, 5 mod 9 and r > 0
        q - q^3 - 2*q^5 + O(q^6) | q - q^3 + 2*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 4
        q - q^3 - 2*q^5 + O(q^6) | q + q^3 + 2*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 4
   The condition that psi == 12 mod 16 and the condition that psi == 3, 6 mod 9 and r > 0
        q - q^3 - 2*q^5 + O(q^6) | q - q^3 - 2*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 8
   The condition that psi == 12 mod 16 and the condition that psi == 2, 5 mod 9 and r > 0
        q - q^3 - 2*q^5 + O(q^6) | q - q^3 - 2*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 4
   The condition that psi == 2 mod 4 and the condition that psi == 3, 6 mod 9 and r > 0
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
   The condition that psi == 2 mod 4 and the condition that psi == 2, 5 mod 9 and r > 0
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q - q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q - q^3 + 2*q^5 + O(q^6) | q + q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q - q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 - 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 + O(q^6) | q - 2*q^5 + O(q^6) | 2
        q + q^3 + 2*q^5 + O(q^6) | q + q^3 + 4*q^5 + O(q^6) | q - 2*q^5 + O(q^6) | 2

The last example of Section 6 works similar to the first.

::

   sage: N1 = (13*L1).value()
   sage: N2 = (13*L2).value()
   sage: cases = conditional_product(F, G, N1, N2)
   sage: pAdics = pAdicBase(QQ, 2)
   sage: cases = [(val, TreeCondition(con.pAdic_tree(pAdics=pAdics))) for val, con in cases
   ....:          if not con.pAdic_tree(pAdics=pAdics).is_empty()]

Now we do all the cases at once

::

   sage: for (F_, G_, N1, N2), con in cases:
   ....:     nfs = [(f, g) for f in get_newforms(N1) for g in get_newforms(N2)]
   ....:     nfs = eliminate_by_traces((F_, G_), nfs, primes=[3, 5, 7, 11])
   ....:     if len(nfs) > 0:
   ....:         print(con)
   ....:         print("    ", lcm([n for f, g, n in nfs]))         
   The condition that psi == 0 mod 128
        21
   The condition that psi == 4 mod 16
        4
   The condition that psi == 8 mod 16
        8
   The condition that psi == 12 mod 16
        20
   The condition that psi == 2 mod 4
        120
