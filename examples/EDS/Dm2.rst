Auto generated from Dm2.sage
============================

.. linkall

Setup
-----

::

   sage: from modular_method.diophantine_equations.conditions import apply_to_conditional_value
   sage: load('frey_curves.sage')

Setting up $E_D(\Q)$
--------------------

::

   sage: D = -2
   sage: Em2 = EllipticCurve([0, 0, 0, D, 0])
   sage: P, = Em2.gens(); T = Em2.torsion_points()[0]
   sage: trace_primes = [p for p in prime_range(50) if not p.divides(2*D)]

Determining the points on $E_D(\Q)$ with B not divisible by a prime >3
----------------------------------------------------------------------

::

   sage: S_integral_points = [Em2(P.Coordinates().sage()) for P in magma(Em2).SIntegralPoints([2, 3])]
   sage: S_integral_points += [-P for P in S_integral_points]
   sage: compare = [
   ....:     P,
   ....:     2*P,
   ....:     T,
   ....:     P + T,
   ....:     2*P + T,
   ....:     3*P + T,
   sage: ]
   sage: compare += [-P for P in compare]
   sage: assert set(S_integral_points) == set(compare)

The case a = 1
--------------

Corresponding to points m*P with m even
::

   sage: a = 1
   sage: E1zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

Conductor computation for a = 1
-------------------------------

::

   sage: e2 = apply_to_conditional_value(lambda E: E.conductor_exponent(2, verbose=True), E1zw)
   Taking into account that B is at least a square
   sage: z, w = E1zw[0][0].parameters()
   sage: con_extra = (CongruenceCondition(w^2 - a*z^4, 2^9) |
   ....:              ~CongruenceCondition(w^2 - a*z^4, 2^2))
   sage: e2 = ConditionalValue([
   ....:     (e, con) for e, con in e2
   ....:     if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
   sage: ]); e2
   1 if ('z', 'w') is 1 of 256 possibilities mod 256
   sage: assert all(
   ....:     mod(ww^2 - a*zz^4, 2^9) == 0
   ....:     for zz, ww in e2[0][1].pAdic_tree().give_as_congruence_condition()[0]
   sage: )

Newform computations for a = 1
------------------------------

::

   sage: Enfs1 = apply_to_conditional_value(
   ....:     lambda E: apply_to_conditional_value(
   ....:         lambda nfs: (E, nfs),
   ....:         E.newform_candidates(
   ....:             bad_primes=(2*D).prime_factors(),
   ....:             algorithm='magma',
   ....:         ),
   ....:     ),
   ....:     E1zw,
   sage: )
   sage: Enfs1 = ConditionalValue([
   ....:     (Enfs, con) for Enfs, con in Enfs1
   ....:     if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
   sage: ])
   sage: assert apply_to_conditional_value(lambda Enfs: len(Enfs[1]), Enfs1) == 0

The case a = -1
---------------

Corresponding to points m*P with m odd
::

   sage: a = -1
   sage: Em1zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

Q-curve data computations for a = -1
------------------------------------

::

   sage: Em1zwg = Em1zw.decomposable_twist()
   sage: K = Em1zwg.definition_field()
   sage: assert K == Em1zwg.decomposition_field()

Newform computation for a = -1
------------------------------

::

   sage: nfsm1 = Em1zwg.newform_candidates(bad_primes=K.primes_above(2*D), algorithm='magma')
   sage: assert apply_to_conditional_value(len, nfsm1) == 16
   sage: z, w = Em1zwg.parameters()
   sage: z, w = z.change_ring(QQ), w.change_ring(QQ)
   sage: nfsm1 = eliminate_by_traces(
   ....:     Em1zwg,
   ....:     nfsm1,
   ....:     condition=CoprimeCondition([z, w]),
   ....:     primes=trace_primes,
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: sum(1 for nf in nfs if nf[-1] == 0),
   ....:     nfsm1,
   sage: ) == 4
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(),
   ....:     nfsm1,
   sage: ) == [2, 3, 7]

Considering odd multiples of P1 = 3*P
-------------------------------------

::

   sage: P1 = 3*P; P1.xy()
   (-1/169, 239/2197)
   sage: assert P1.xy()[0].denominator().prime_factors() == [13]
   sage: nfsm1P = eliminate_by_trace(
   ....:     Em1zwg,
   ....:     nfsm1,
   ....:     13,
   ....:     condition=(CoprimeCondition([z, w]) &
   ....:                CongruenceCondition(w^2 - a*z^4, 13)),
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
   ....:     nfsm1P
   sage: ) == [2, 5, 7, 13]
   l-th power must be multiple of 13*P1
   sage: assert [
   ....:     p for p in prime_range(100)
   ....:     if p.divides((13*P1).xy()[0].denominator())
   sage: ] == [13, 37, 41]
   sage: nfsm1P = eliminate_by_trace(
   ....:     Em1zwg,
   ....:     nfsm1P,
   ....:     37,
   ....:     condition=(CoprimeCondition([z, w]) &
   ....:                CongruenceCondition(w^2 - a*z^4, 37)),
   ....:     verbose=True,
   sage: )
   sage: nfsm1P = eliminate_by_trace(
   ....:     Em1zwg,
   ....:     nfsm1P,
   ....:     41,
   ....:     condition=(CoprimeCondition([z, w]) &
   ....:                CongruenceCondition(w^2 - a*z^4, 41)),
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
   ....:     nfsm1P
   sage: ) == [2]

The case a = 2
--------------

Corresponding to points m*P + T with m odd
::

   sage: a = 2
   sage: E2zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

Q-curve data computations for a = 2
-----------------------------------

::

   sage: E2zwg = E2zw.decomposable_twist()
   sage: K = E2zwg.definition_field()
   sage: assert K == E2zwg.decomposition_field()

Newform computation for a = 2
-----------------------------

::

   sage: nfs2 = E2zwg.newform_candidates(bad_primes=K.primes_above(2*D), algorithm='magma')
   sage: assert apply_to_conditional_value(len, nfs2) == 28
   sage: z, w = E2zwg.parameters()
   sage: z, w = z.change_ring(QQ), w.change_ring(QQ)
   sage: nfs2 = eliminate_by_traces(
   ....:     E2zwg,
   ....:     nfs2,
   ....:     condition=CoprimeCondition([z, w]),
   ....:     primes=trace_primes,
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: sum(1 for nf in nfs if nf[-1] == 0),
   ....:     nfs2,
   sage: ) == 12
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(),
   ....:     nfs2,
   sage: ) == [2]

Considering odd multiples of P1 = 5*P + T
-----------------------------------------

::

   sage: P1 = 5*P + T; P1.xy()
   (4651250/1803649, -8388283850/2422300607)
   sage: assert P1.xy()[0].denominator().prime_factors() == [17, 79]
   sage: nfs2P = eliminate_by_trace(
   ....:     E2zwg,
   ....:     nfs2,
   ....:     17,
   ....:     condition=(CoprimeCondition([z, w]) &
   ....:                CongruenceCondition(w^2 - a*z^4, 17)),
   ....:     verbose=True,
   sage: )
   sage: nfs2P = eliminate_by_trace(
   ....:     E2zwg,
   ....:     nfs2,
   ....:     79,
   ....:     condition=(CoprimeCondition([z, w]) &
   ....:                CongruenceCondition(w^2 - a*z^4, 79)),
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
   ....:     nfs2P
   sage: ) == [2, 5, 79]

The case a = -2
---------------

Corresponding to points m*P + T with m even
::

   sage: a = -2
   sage: Em2zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

Q-curve data computations for a = -2
------------------------------------

::

   sage: Em2zwg = Em2zw.decomposable_twist()
   sage: K = Em2zwg.definition_field()
   sage: assert K == Em2zwg.decomposition_field()

Newform computation for a = -2
------------------------------

::

   sage: nfsm2 = Em2zwg.newform_candidates(bad_primes=K.primes_above(2*D), algorithm='magma')
   sage: assert apply_to_conditional_value(len, nfsm2) == 28
   sage: z, w = Em2zwg.parameters()
   sage: z, w = z.change_ring(QQ), w.change_ring(QQ)
   sage: nfsm2 = eliminate_by_traces(
   ....:     Em2zwg,
   ....:     nfsm2,
   ....:     condition=CoprimeCondition([z, w]),
   ....:     primes=trace_primes,
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: sum(1 for nf in nfs if nf[-1] == 0),
   ....:     nfsm2,
   sage: ) == 4
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(),
   ....:     nfsm2,
   sage: ) == [2, 3, 7]

Considering odd multiples of P1 = 2*P + T
-----------------------------------------

::

   sage: P1 = 2*P + T; P1.xy()
   (-8/9, -28/27)
   sage: assert P1.xy()[0].denominator().prime_factors() == [3]
   sage: nfsm2P = eliminate_by_trace(
   ....:     Em2zwg,
   ....:     nfsm2,
   ....:     3,
   ....:     condition=(CoprimeCondition([z, w]) &
   ....:                CongruenceCondition(w^2 - a*z^4, 3)),
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
   ....:     nfsm2P
   sage: ) == [2, 3, 7]
   l-th power must be multiple of 3*P1
   sage: assert [
   ....:     p for p in prime_range(100)
   ....:     if p.divides((3*P1).xy()[0].denominator())
   sage: ] == [3, 11]
   sage: nfsm2P = eliminate_by_trace(
   ....:     Em2zwg,
   ....:     nfsm2P,
   ....:     11,
   ....:     condition=(CoprimeCondition([z, w]) &
   ....:                CongruenceCondition(w^2 - a*z^4, 11)),
   ....:     verbose=True,
   sage: )
   sage: assert apply_to_conditional_value(
   ....:     lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
   ....:     nfsm2P
   sage: ) == [2, 3]