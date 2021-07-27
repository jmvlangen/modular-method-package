# This file introduces the function
# Frey_curve_of_divisibility_sequence that allows to compute the curve
# E_{a, z, w} from the values of a and D The returned curve will be
# FreyCurve if a is a square and a FreyQcurve otherwise.  Note that
# the curve is initialized with the condition that (z, w) are coprime
# and that there exists an integer B such that a*z^4 + (D/a)*B^4 ==
# w^2. Note that the latter is only enforced up to a given precision
# (default 1, as computation time drastically increases
# otherwise). For the prime above 2 a seperate condition is added that
# is the previously mentioned condition up to a precision 3 higher
# than the given precision.
#
# Whenever a == 1 the function Frey_curve_of_divisibility_sequence_1
# is used, which precomputes the conductor exponent at 2 using the
# results from Table 5.1. Furthermore its result is a conditional
# value containing the curve with the lowest possible conductor
# exponent at 2 for each case.

from modular_method import *
from modular_method.number_fields.field_constructors import field_with_root
from modular_method.diophantine_equations.conditions import ConditionalValue
from modular_method.diophantine_equations.conditions import TreeCondition
from modular_method.padics.pAdic_base import pAdicBase

def Frey_curve_of_divisibility_sequence_1(R, C, C2):
    z, w = R.gens()
    pAdics = pAdicBase(QQ, 2)
    conditions = {1: {}, -1: {}, 2: {}, -2: {}}
    # conditions[gamma][e] will be the condition to get
    # conductor exponent e with the twist gamma
    # All conditions are chosen to be exclusive
    conditions[1][8] = TreeCondition((
        C2 &
        ~CongruenceCondition(w^2 - z^4, 2)
    ).pAdic_tree(pAdics=pAdics))
    conditions[1][7] = TreeCondition((
        C2 &
        CongruenceCondition(w^2 - z^4 - 8, 16)
    ).pAdic_tree(pAdics=pAdics))
    conditions[-1][3] = TreeCondition((
        C2 &
        ((CongruenceCondition(z - 1, 4) &
          CongruenceCondition(w + z^2 - 8, 32)) |
         (CongruenceCondition(z - 3, 4) &
          CongruenceCondition(w + z^2, 2^5) &
          ~CongruenceCondition(w + z^2, 2^7)))
    ).pAdic_tree(pAdics=pAdics))
    conditions[1][3] = TreeCondition((
        C2 &
        ((CongruenceCondition(z - 3, 4) &
          CongruenceCondition(w + z^2 - 8, 32)) |
         (CongruenceCondition(z - 1, 4) &
          CongruenceCondition(w + z^2, 2^5) &
          ~CongruenceCondition(w + z^2, 2^7)))
    ).pAdic_tree(pAdics=pAdics))
    conditions[1][2] = TreeCondition((
        C2 &
        CongruenceCondition(z - 1, 4) &
        CongruenceCondition(w + z^2 - 24, 32)
    ).pAdic_tree(pAdics=pAdics))
    conditions[-1][2] = TreeCondition((
        C2 &
        CongruenceCondition(z - 3, 4) &
        CongruenceCondition(w + z^2 - 24, 32)
    ).pAdic_tree(pAdics=pAdics))
    conditions[-2][2] = TreeCondition((
        C2 &
        CongruenceCondition(z - 1, 4) &
        CongruenceCondition(w - z^2 - 8, 32)
    ).pAdic_tree(pAdics=pAdics))
    conditions[2][2] = TreeCondition((
        C2 &
        CongruenceCondition(z - 3, 4) &
        CongruenceCondition(w - z^2 - 8, 32)
    ).pAdic_tree(pAdics=pAdics))
    conditions[2][3] = TreeCondition((
        C2 &
        ((CongruenceCondition(z - 1, 4) &
          CongruenceCondition(w - z^2 - 24, 32)) |
         (CongruenceCondition(z - 3, 4) &
          CongruenceCondition(w - z^2, 2^5) &
          ~CongruenceCondition(w - z^2, 2^7)))
    ).pAdic_tree(pAdics=pAdics))
    conditions[-2][3] = TreeCondition((
        C2 &
        ((CongruenceCondition(z - 3, 4) &
          CongruenceCondition(w - z^2 - 24, 32)) |
         (CongruenceCondition(z - 1, 4) &
          CongruenceCondition(w - z^2, 2^5) &
          ~CongruenceCondition(w - z^2, 2^7)))
    ).pAdic_tree(pAdics=pAdics))
    conditions[1][5] = TreeCondition((
        C2 &
        CongruenceCondition(w + z^2 - 16, 32)
    ).pAdic_tree(pAdics=pAdics))
    conditions[2][5] = TreeCondition((
        C2 &
        CongruenceCondition(w - z^2 - 16, 32)
    ).pAdic_tree(pAdics=pAdics))
    conditions[1][0] = TreeCondition((
        C2 &
        CongruenceCondition(z - 1, 4) &
        CongruenceCondition(w + z^2 - 2^7, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    conditions[-1][0] = TreeCondition((
        C2 &
        CongruenceCondition(z - 3, 4) &
        CongruenceCondition(w + z^2 - 2^7, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    conditions[-2][0] = TreeCondition((
        C2 &
        CongruenceCondition(z - 1, 4) &
        CongruenceCondition(w - z^2 - 2^7, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    conditions[2][0] = TreeCondition((
        C2 &
        CongruenceCondition(z - 3, 4) &
        CongruenceCondition(w - z^2 - 2^7, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    conditions[1][1] = TreeCondition((
        C2 &
        CongruenceCondition(z - 1, 4) &
        CongruenceCondition(w + z^2, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    conditions[-1][1] = TreeCondition((
        C2 &
        CongruenceCondition(z - 3, 4) &
        CongruenceCondition(w + z^2, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    conditions[-2][1] = TreeCondition((
        C2 &
        CongruenceCondition(z - 1, 4) &
        CongruenceCondition(w - z^2, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    conditions[2][1] = TreeCondition((
        C2 &
        CongruenceCondition(z - 3, 4) &
        CongruenceCondition(w - z^2, 2^8)
    ).pAdic_tree(pAdics=pAdics))
    result = []
    for gamma in conditions:
        conE2 = reduce(lambda C1, C2: C1 | C2, conditions[gamma].values())
        conE = C & conE2
        E = FreyCurve(
            [0, 4*z*gamma, 0, 2*(z^2 + w)*gamma^2, 0],
            condition=conE
        )
        N = [(val, con) for val, con in conditions[gamma].items() if not con.never()]
        if len(N) > 0:
            if len(N) == 1:
                E.conductor_exponent.set_cache(N[0], 2)
            else:
                E.conductor_exponent.set_cache(ConditionalValue(N), 2)
            result.append((E, conE2))
    if len(result) > 1:
        return ConditionalValue(result)
    else:
        return result[0]

def Frey_curve_of_divisibility_sequence(a, D, precision=1):
    R.<z, w, B> = QQ[]
    d = ZZ(D / a)
    C = ExistsCondition(a*z^4 + d*B^4 - w^2, [z, w], precision=precision)
    C2 = TreeCondition((CoprimeCondition([z, w]) &
                        ExistsCondition(
                            a*z^4 + d*B^4 - w^2,
                            [z, w],
                            precision=precision+3,
                        )).pAdic_tree(pAdics=pAdicBase(QQ, 2)))
    z, w = QQ[z, w].gens()
    C = CoprimeCondition([z, w]) & C
    if a == 1:
        return Frey_curve_of_divisibility_sequence_1(z.parent(), C, C2)
    C = C & C2
    K = field_with_root(QQ, a)
    sqrta = sqrt(K(a))
    invariants = [0, 4*sqrta*z, 0, 2*(a * z^2 + sqrta * w), 0]
    if K == QQ:
        return FreyCurve(invariants, condition=C)
    return FreyQcurve(invariants, condition=C, guessed_degrees=[2])
