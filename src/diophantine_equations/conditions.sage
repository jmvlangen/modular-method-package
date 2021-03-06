r"""Conditions on variables and functions to work with them.

Gives some standard diophantine conditions that one might want to
impose on a collection of variables/parameters.

Simple examples of conditions on variables are polynomial conditions
and coprimality conditions. Other implemented conditions include the
condition for a polynomial to be zero modulo an ideal and the
condition for a polynomial to be an n-th power. Conditions can also be
combined, by taking their negation or a binary AND or OR.

Each condition can also be turned into a pAdicTree which collects all
the possible p-adic values for the variables for which the condition
might hold. Each tree also again defines a condition known as a
TreeCondition.

Using the class ConditionalValue one can work with values which might
depend on some condition to hold. These values can be used in
arithmetic expressions which can again be evaluated to conditional
values themselves. One can also use ConditionalValues in functions by
using the method :func:`apply_to_conditional_value`.

EXAMPLES::

    sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
    sage: R.<x,y> = QQ[]
    sage: CoprimeCondition([x, y])
    The condition that the variables ('x', 'y') are pairwise coprime.
    sage: PolynomialCondition(x^2 + y^2 - 4)
    The condition that x^2 + y^2 - 4 == 0

AUTHORS:

- Joey van Langen (2019-02-13): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from functools import reduce

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.functions.other import floor, ceil

from sage.all import Integer, RealNumber, QQ, ZZ

from sage.rings.polynomial.multi_polynomial import is_MPolynomial
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.rings.ideal import is_Ideal

from modular_method.padics.pAdic_tree import pAdicTree
from modular_method.padics.pAdic_solver import find_pAdic_roots

class Condition_base(SageObject):
    """A class that represents a condition on some variables.

    This class should not be used by itself.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
        sage: R.<x,y> = QQ[]
        sage: CoprimeCondition([x, y])
        The condition that the variables ('x', 'y') are pairwise coprime.
        sage: PolynomialCondition(x^2 + y^2 - 4)
        The condition that x^2 + y^2 - 4 == 0

    """

    def __init__(self, variables):
        r"""Initialize a Condition.

        INPUT:
        
        - ``variables`` -- A list or tuple of variables on which this
          condition applies. Each entry of the given list or tuple
          will be converted to a string.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<x,y> = QQ[]
            sage: CoprimeCondition([x, y])
            The condition that the variables ('x', 'y') are pairwise coprime.
            sage: CoprimeCondition(['x', 'y'])
            The condition that the variables ('x', 'y') are pairwise coprime.

        """
        self._vars = tuple(str(v) for v in variables)

    def variables(self):
        r"""Give the variables on which this Condition applies

        OUTPUT:

        A tuple of variables on which this Condition applies. Each
        variable is represented as a string.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<x, y> = QQ[]
            sage: CoprimeCondition([x, y]).variables()
            ('x', 'y')

        Note that even when a ring has multiple variables only the
        relevant variables for this condition are returned::

            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y, z> = ZZ[]
            sage: PolynomialCondition(x^2 - y^2 - 4).variables()
            ('x', 'y')

        """
        return self._vars

    def never(self):
        r"""Tell if this condition never holds

        .. NOTE::

        This function returning False does not imply this condition
        may not hold, as complex conditions might not be able to
        determine whether they never hold or not.

        OUTPUT:

        True or False. If the return value is True this condition can
        not hold on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.never()
            True
            sage: C2.never()
            False

        Note that when combining conditions, a condition that never
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that never holds
            sage: C1 | C2
            The condition that the variables ('a', 'b') are pairwise coprime.

        This method is not able to know of every combination of
        conditions whether it never holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 & C2).never()
            False

        """
        return False

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        return False
        
    def __and__(self, other):
        r"""Create the condition that both conditions hold.

        INPUT:
        
        - ``other`` -- A Condition, i.e. an instance of
          Condition_base.

        OUTPUT:

        A Condition object that holds on all values where both this
        Condition object and the given Condition object hold.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: CoprimeCondition([x,y]) & PolynomialCondition(x^2 + y^2 - 4)
            The condition that the variables ('x', 'y') are pairwise coprime and the condition that x^2 + y^2 - 4 == 0

        If two conditions are the same, the 'and' of both of them is
        just the first::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = ZZ[]
            sage: C1 = CongruenceCondition(x, 3); C1
            The condition that x == 0 modulo 3
            sage: C2 = CongruenceCondition(x, 3); C2
            The condition that x == 0 modulo 3
            sage: C1 & C2
            The condition that x == 0 modulo 3

        .. SEEALSO ::
        
           :class:`AndCondition`

        """
        if (self == other or self.never() or other == None or
            other.always()):
            return self
        elif other.never() or self.always():
            return other
        else:
            return AndCondition(self, other)

    def __or__(self, other):
        r"""Create the condition that either condition holds.

        INPUT:
        
        - ``other`` -- A Condition, i.e. an instance of
          Condition_base.

        OUTPUT:

        A Condition object that holds on all values where either this
        Condition object or the given Condition object holds.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: CoprimeCondition([x,y]) | PolynomialCondition(x^2 + y^2 - 4)
            The condition that the variables ('x', 'y') are pairwise coprime or the condition that x^2 + y^2 - 4 == 0

        If two conditions are the same, the 'or' of both of them is
        just the first::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = ZZ[]
            sage: C1 = CongruenceCondition(x, 3); C1
            The condition that x == 0 modulo 3
            sage: C2 = CongruenceCondition(x, 3); C2
            The condition that x == 0 modulo 3
            sage: C1 | C2
            The condition that x == 0 modulo 3

        .. SEEALSO::

            :class:`OrCondition`

        """
        if (self == other or self.always() or other == None or
            other.never()):
            return self
        elif self.never() or other.always():
            return other
        else:
            return OrCondition(self, other)

    def __invert__(self):
        r"""Create the condition that this condition does not hold.

        OUTPUT:

        A Condition object that holds on all values where this
        Condition does not hold.

        EXAMPLES::
        
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: ~ PolynomialCondition(x^2 + y^2 - 4)
            The condition that x^2 + y^2 - 4 ~= 0

        Note that a double not simplifies in print, but gives a
        different object::

            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1); C
            The condition that -x^3 + y^2 - 1 == 0
            sage: ~~C
            The condition that -x^3 + y^2 - 1 == 0
            sage: C == ~~C
            False

        .. SEEALSO::

            :class:`NotCondition`

        """
        return NotCondition(self)

    def _repr_(self):
        return "A condition on the variables %s"%(self._vars,)

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree that satisfies the
        condition defined by this object on the variables therein.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree which could satisfy the condition defined by this
        object. If complement was set to True will return a tuple with
        the afore mentioned as its first entry and a pAdicTree
        containing those values of the given pAdicTree for which this
        condition could not be satisfied.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        raise NotImplementedError("The method 'pAdic_tree' is not " +
                                  "implemented for this condition class")

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                isinstance(self, other.__class__) and
                self.variables() == other.variables())

    def __ne__(self, other):
        return not self.__eq__(other)

class PolynomialCondition(Condition_base):
    r"""The condition that a certain polynomial is zero.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
        sage: R.<x, y> = ZZ[]
        sage: PolynomialCondition(x^2 + y^2 - 4)
        The condition that x^2 + y^2 - 4 == 0
        sage: PolynomialCondition(y^2 - x^3 - 1)
        The condition that -x^3 + y^2 - 1 == 0

    """

    def __init__(self, polynomial, precision=20):
        r"""Initialize a PolynomialCondition.

        INPUT:

        - ``polynomial`` -- A polynomial in one or more variables
          which should be zero

        - ``precision`` -- A strictly positive integer (default 20)
          which determines the default p-adic precision used when
          computing p-adic trees.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: PolynomialCondition(x^2 + y^2 - 4)
            The condition that x^2 + y^2 - 4 == 0
            sage: PolynomialCondition(y^2 - x^3 - 1)
            The condition that -x^3 + y^2 - 1 == 0

        """
        if not (is_Polynomial(polynomial) or is_MPolynomial(polynomial)):
            raise ValueError("The given argument " + str(polynomial) +
                             " is not a polynomial.")
        self._f = polynomial
        self._prec = precision
        Condition_base.__init__(self, self._f.variables())

    def never(self):
        r"""Tell if this condition never holds

        .. NOTE::

        This function returning False does not imply this condition
        may not hold, as complex conditions might not be able to
        determine whether they never hold or not.

        OUTPUT:

        True or False. If the return value is True this condition can
        not hold on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.never()
            True
            sage: C2.never()
            False

        Note that when combining conditions, a condition that never
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that never holds
            sage: C1 | C2
            The condition that the variables ('a', 'b') are pairwise coprime.

        This method is not able to know of every combination of
        conditions whether it never holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 & C2).never()
            False

        """
        return (self._f.is_constant() and self._f != 0)

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        return self._f == 0

    def polynomial(self):
        r"""Give the polynomial associated to this condition.

        OUTPUT:

        The polynomial that defines this PolynomialCondition.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: PolynomialCondition(x^2 + y^2 - 4).polynomial()
            x^2 + y^2 - 4

        """
        return self._f

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   precision=None, verbose=False, **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those values for
        the variables such that the polynomial of this condition is
        zero on them.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        - ``precision`` -- A strictly positive integer giving up to
          what precision the resulting tree should be found. If
          unspecified will be set to the precision given upon
          initialization of this object.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings. A larger value will lead to
          more information being printed.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree which satisfies the polynomial of this condition
        being equal to zero modulo P^n, where P is the prime defined
        by the given pAdics and n is the given precision. If
        complement was set to True will return a tuple with the afore
        mentioned as its first entry and the complement of that tree
        within the given pAdicTree as its second argument.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if precision == None:
            precision = self._prec
        if "precision_cap" in kwds and precision > kwds["precision_cap"]:
            if verbose >= 0:
                print("Warning: A p-Adic tree for " + str(self) +
                       " is computed with a lower precision than the given " +
                       "precision " + str(precision) + ", due to a given " +
                       "precision_cap " + str(kwds["precision_cap"]) + "!")
            precision = kwds["precision_cap"]
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        big_vars = list(pAdic_tree.variables())
        for var in self.variables():
            if var not in big_vars:
                big_vars.append(var)
        K = pAdics.number_field()
        K0 = self.polynomial().parent().base()
        if K0.is_subring(QQ):
            iota = K0.hom(K)
        else:
            iota = K0.hom([a.minpoly().change_ring(K).roots()[0][0] for a in K0.gens()], K)
        S = PolynomialRing(K, big_vars)
        Tyes, Tno = find_pAdic_roots(S(self.polynomial().change_ring(iota)),
                                    pAdics=pAdics,
                                    variables=[S(v) for v in pAdic_tree.variables()],
                                    value_tree=pAdic_tree.root(),
                                    precision=precision,
                                    verbose=verbose)
        Tyes = pAdicTree(variables=big_vars, root=Tyes)
        Tyes = Tyes.change_variables_to(pAdic_tree.variables())
        if complement:
            Tno = pAdicTree(variables=big_vars, root=Tno)
            Tno = Tno.change_variables_to(pAdic_tree.variables())
            return Tyes, Tno
        else:
            return Tyes

    def _repr_(self):
        return "The condition that %s == 0"%(self.polynomial(),)

    def _repr_short(self):
        return "%s == 0"%(self.polynomial(),)
        
    def _latex_(self):
        return latex(self.polynomial()) + " = " + latex(0)

    def _cache_key(self):
        return 'TreeCondition', self.polynomial()

    def __eq__(self, other):
        return (Condition_base.__eq__(self, other) and
                self.polynomial() == other.polynomial())

class ExistsCondition(PolynomialCondition):
    r"""The condition that some additional variables satisfying a
    polynomial relation exist

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import ExistsCondition
        sage: R.<x, y> = ZZ[]
        sage: ExistsCondition(x^3 - y^2, [x])
        The condition that x^3 - y^2 == 0 for some y

    An ExistsCondition could be used to assert that some variable is a
    square::

        sage: from modular_method.padics.pAdic_base import pAdicBase
        sage: from modular_method.diophantine_equations.conditions import ExistsCondition
        sage: from modular_method.diophantine_equations.conditions import TreeCondition
        sage: R.<x, y> = QQ[]
        sage: C = ExistsCondition(x - y^2, [x], precision=4); C
        The condition that -y^2 + x == 0 for some y
        sage: TreeCondition(C.pAdic_tree(pAdics=pAdicBase(QQ, 2)))
        The condition that x == 0, 1, 4, 9 mod 16
        sage: TreeCondition(C.pAdic_tree(pAdics=pAdicBase(QQ, 3)))
        The condition that x == 0, 1, 4, 7, 9, 10, 13, 16, 19, 22, 25, 28, 31, 34, 36, 37, 40, 43, 46, 49, 52, 55, 58, 61, 63, 64, 67, 70, 73, 76, 79 mod 81

    """

    def __init__(self, polynomial, variables, precision=20):
        r"""Initialize a ExistsCondition

        INPUT:

        - ``polynomial`` -- A polynomial that should be satisfied

        - ``variables`` -- The variables of this condition. All other
          variables will be considered as bound variables bound by an
          existential quantifier.

        - ``precision`` -- A strictly positive integer (default: 20)
          which will be the default p-adic precision used when
          computing p-adic trees for this condition.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import ExistsCondition
            sage: R.<x, y, z> = ZZ[]
            sage: ExistsCondition(x^2 + y^2 + z^2, [x])
            The condition that x^2 + y^2 + z^2 == 0 for some y, and z
            sage: ExistsCondition(x^2 + y^2 + z^2, [x, y])
            The condition that x^2 + y^2 + z^2 == 0 for some z

        """
        PolynomialCondition.__init__(self, polynomial, precision=precision)
        allvars = list(self._vars)
        self._vars = tuple(str(var) for var in variables)
        for var in self._vars:
            try:
                allvars.remove(var)
            except ValueError:
                raise ValueError("Variable '" + var + "' is not among " +
                                 "the variables in the polynomial or is " +
                                 "a duplicate")
        self._bound_vars = tuple(allvars)

    def bound_variables(self):
        r"""Give the bound variables of this condition

        OUTPUT:

        A tuple of bound variables in this Condition. Each variable is
        represented as a string.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import ExistsCondition
            sage: R.<a, b, c> = ZZ[]
            sage: C = ExistsCondition((a - b)^2 + c^3, [a, b])
            sage: C.bound_variables()
            ('c',)
            sage: C = ExistsCondition((a - b)^2 + c^3, [c])
            sage: C.bound_variables()
            ('a', 'b')

        Note that only variables appearing in the polynomial of this
        condition can be bound::

            sage: from modular_method.diophantine_equations.conditions import ExistsCondition
            sage: R.<x, y, z> = ZZ[]
            sage: C = ExistsCondition(x^2 - y^3, [x])
            sage: C.bound_variables()
            ('y',)

        .. SEE_ALSO::
        
            :meth:`variables`
            :meth:`polynomial`

        """
        return self._bound_vars

    def pAdic_tree(self, pAdic_tree=None, pAdics=None,
                   complement=False, precision=None, verbose=False,
                   **kwds):
        r"""Give this condition as a pAdicTree

        Given a p-adic tree, returns the subtree of those values for
        the variables for which there exist values for the bound
        variables such that the polynomial of this condition is
        satisfied.

        INPUT:

        - ``pAdic_tree`` -- A pAdicTree object (default: None) on
          which this condition should be applied. If set to None will
          be initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should also be returned.

        - ``precision`` -- A strictly positive integer giving up to
          what precision the resulting tree should be found. If
          unspecified will be set to the precision given upon
          initialization of this object.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger than zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments. If set to any negative value will also prevent the
          printing of any warnings. A larger value will lead to more
          information being printed.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree for which values of the bound variables exists such
        that the polynomial of this condition is zero modulo P^n. Here
        P is the prime defined by the given pAdics and n is the given
        precision.

        If complement was set to True will return a tuple with the
        afore mentioned as its first entry and the complement of that
        tree within the given pAdicTree as its second argument.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if precision == None:
            precision = self._prec
        if "precision_cap" in kwds and precision > kwds["precision_cap"]:
            if verbose >= 0:
                print("Warning: A p-Adic tree for " + str(self) +
                       " is computed with a lower precision than the given " +
                       "precision " + str(precision) + ", due to a given " +
                       "precision_cap " + str(kwds["precision_cap"]) + "!")
            precision = kwds["precision_cap"]
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        for var in self._vars:
            if not (var in pAdic_tree.variables()):
                raise ValueError("Variable '" + var + "' not part of " +
                                 "the given p-adic tree.")
        for var in self._bound_vars:
            if var in pAdic_tree.variables():
                raise ValueError("Bound variable '" + var + "' part of " +
                                 "the given p-adic tree.")
        big_tree = pAdic_tree.add_variable(*self._bound_vars)
        K = pAdics.number_field()
        K0 = self.polynomial().parent().base()
        if K0.is_subring(QQ):
            iota = K0.hom(K)
        else:
            iota = K0.hom([a.minpoly().change_ring(K).roots()[0][0] for a in K0.gens()], K)
        S = PolynomialRing(K, big_tree.variables())
        Tyes, Tno = find_pAdic_roots(S(self.polynomial().change_ring(iota)),
                                     pAdics=pAdics,
                                     variables=[S(v) for v in big_tree.variables()],
                                     value_tree=big_tree.root(),
                                     precision=precision,
                                     verbose=verbose)
        Tyes = pAdicTree(variables=big_tree.variables(), root=Tyes)
        Tyes = Tyes.change_variables_to(pAdic_tree.variables())
        if complement:
            Tno = pAdic_tree.difference(Tyes)
            return Tyes, Tno
        else:
            return Tyes

    def _repr_(self):
        result = PolynomialCondition._repr_(self)
        n = len(self._bound_vars)
        if n > 0:
            result += " for some "
            for i in range(n):
                if i > 0 and i == n - 1:
                    result += "and "
                result += self._bound_vars[i]
                if i < n - 1:
                    result += ", "
        return result

    def _repr_short(self):
        result = PolynomialCondition._repr_short(self)
        n = len(self._bound_vars)
        if n > 0:
            result += " for some "
            for i in range(n):
                if i > 0 and i == n - 1:
                    result += "and "
                result += self._bound_vars[i]
                if i < n - 1:
                    result += ", "
        return result

    def _latex_(self):
        result = ""
        n = len(self._bound_vars)
        if n > 0:
            result += "\\exists "
            for i in range(n):
                result += self._bound_vars[i]
                if i < n - 1:
                    result += ", "
            result += " : "
        result += PolynomialCondition._latex_(self)
        return result

    def _cache_key(self):
        return 'ExistsCondition', self.polynomial(), self.variables()

    def __eq__(self, other):
        return (Condition_base.__eq__(self, other) and
                self.polynomial() == other.polynomial())

class CongruenceCondition(PolynomialCondition):
    r"""The condition that a polynomial is congruent to zero.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
        sage: R.<x, y> = ZZ[]
        sage: CongruenceCondition(x, 3)
        The condition that x == 0 modulo 3
        sage: CongruenceCondition(y^2 - 4, 12)
        The condition that y^2 - 4 == 0 modulo 12

    """

    def __init__(self, polynomial, modulus):
        r"""Initialize a CongruenceCondition.

        INPUT:

        - ``polynomial`` -- A polynomial for which a congruence should
          hold.

        - ``modulus`` -- An algebraic integer or integral ideal of the
          ring of integers of a number field modulo which the
          polynomial should be zero.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: CongruenceCondition(x, 3)
            The condition that x == 0 modulo 3
            sage: CongruenceCondition(y^2 - 4, 12)
            The condition that y^2 - 4 == 0 modulo 12

        We can also work over some number field::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: K.<w> = QuadraticField(30)
            sage: S.<x, y> = K[]
            sage: I = K.prime_above(3)^5
            sage: CongruenceCondition(x^2 + y^2, I)
            The condition that x^2 + y^2 == 0 modulo (27, 9*w)

        """
        PolynomialCondition.__init__(self, polynomial)
        self._mod = modulus
        if hasattr(self._mod, 'is_principal') and self._mod.is_principal():
            self._mod = self._mod.gens_reduced()[0]    

    def modulus(self):
        r"""Give the modulus considered in this condition.

        OUTPUT:

        The algebraic integer or integral ideal modulo which the
        polynomial of this condition should be zero.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: CongruenceCondition(x^2 + 2*y, 7).modulus()
            7

        Over number fields the modulus can be an ideal::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: K.<w> = QuadraticField(30)
            sage: S.<x, y> = K[]
            sage: I = K.prime_above(3)^5
            sage: CongruenceCondition(x^2 + y^2, I).modulus()
            Fractional ideal (27, 9*w)

        However, if the ideal is principal it will be replaced by a
        generator::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: K.<w> = QuadraticField(30)
            sage: S.<x, y> = K[]
            sage: I = K.prime_above(3)^4
            sage: CongruenceCondition(x^2 + y^2, I).modulus()
            9

        """
        return self._mod
        
    def pAdic_tree(self, pAdic_tree=None, pAdics=None,
                   complement=False, verbose=False, **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those values for
        the variables such that the polynomial of this condition can
        be zero modulo the modulus given in this condition.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings. A larger value will lead to
          more information being printed.

        OUTPUT:

        If the prime given by the given pAdics divides the modulus of
        this Condition, returns the pAdicTree containing all values
        for the variables in the give pAdicTree such that the
        polynomial of this condition is zero modulo the maximal power
        to which this prime ideal appears in the modulus.

        If the prime given does not divide the modulus of this
        Condition the returned tree is simply the one given.

        If complement was set to True will return a tuple with the
        afore mentioned as its first entry. The second entry will be
        the complement of the first pAdicTree inside the given
        pAdicTree if the prime given divides the modulus and the given
        pAdicTree otherwise.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        precision = pAdics.valuation(self.modulus())
        result = PolynomialCondition.pAdic_tree(self, pAdic_tree=pAdic_tree,
                                                pAdics=pAdics,
                                                complement=complement,
                                                verbose=verbose,
                                                precision=precision,
                                                **kwds)
        if complement and precision == 0:
            return result[0], result[0]
        else:
            return result

    def never(self):
        r"""Tell if this condition never holds

        .. NOTE::

        This function returning False does not imply this condition
        may not hold, as complex conditions might not be able to
        determine whether they never hold or not.

        OUTPUT:

        True or False. If the return value is True this condition can
        not hold on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.never()
            True
            sage: C2.never()
            False

        Note that when combining conditions, a condition that never
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that never holds
            sage: C1 | C2
            The condition that the variables ('a', 'b') are pairwise coprime.

        This method is not able to know of every combination of
        conditions whether it never holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 & C2).never()
            False

        TESTS:

        Check that the bug where congruence conditions with
        polynomials without a non-constant term were both always and
        never is fixed::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C = CongruenceCondition(3*x, 3)
            sage: C.always()
            True
            sage: C.never()
            False

        """
        m1 = self._f.parent()(1)
        return (not self._mod.divides(self._f.monomial_coefficient(m1)) and
                all(self._mod.divides(self._f.monomial_coefficient(m))
                    for m in self._f.monomials() if m != m1))

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        if is_Ideal(self._mod):
            return all(cf in self._mod for cf in self._f.coefficients())
        else:
            return all((cf / self._mod).is_integral()
                       for cf in self._f.coefficients())
    
    def _repr_(self):
        mod = self.modulus()
        mod_str = (mod._repr_short() if hasattr(mod, '_repr_short')
                   else str(mod))
        return "The condition that %s == 0 modulo %s"%(self.polynomial(),
                                                       mod_str)

    def _repr_short(self):
        mod = self.modulus()
        mod_str = (mod._repr_short() if hasattr(mod, '_repr_short')
                   else str(mod))
        return "%s == 0 mod %s"%(self.polynomial(),
                                 mod_str)
        
    def _latex_(self):
        return (latex(self.polynomial()) + " = " + latex(0) +
                "\\text{ (mod }" + latex(self.modulus()) + "\\text{)}")

    def _cache_key(self):
        return 'CongruenceCondition', self.polynomial(), self.modulus()

    def __eq__(self, other):
        return (PolynomialCondition.__eq__(self, other) and
                self.modulus() == other.modulus())

class PowerCondition(PolynomialCondition):
    r"""The condition that a polynomial is an n-th power.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import PowerCondition
        sage: R.<x, y> = ZZ[]
        sage: PowerCondition(x^3 + y^3, 3)
        The condition that x^3 + y^3 == x0^n0 with n0 >= 3

    """

    def __init__(self, polynomial, least_exp=1):
        r"""Initialize a PowerCondition.

        INPUT:

        - ``polynomial`` -- A polynomial that should be an unknown
          power of some number.

        - ``least_exp`` -- A strictly positive integer (default: 1)
          that is the least power that this polynomial must be.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import PowerCondition
            sage: R.<x, y> = ZZ[]
            sage: PowerCondition(x^3 + y^3, 3)
            The condition that x^3 + y^3 == x0^n0 with n0 >= 3

        """
        PolynomialCondition.__init__(self, polynomial)
        self._exp = least_exp

    def least_exponent(self):
        r"""Give the least n such that the polynomial is an n-th power.

        OUTPUT:

        The smallest integer such that the polynomial in this
        condition is allowed to be that power of some number.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import PowerCondition
            sage: R.<x, y> = ZZ[]
            sage: PowerCondition(x^3 + y^3, 3).least_exponent()
            3

        """
        return self._exp

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   verbose=False, **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those values for
        the variables such that the polynomial of this condition could
        be at least the power n of some number, where n is the least
        exponent stored in this condition.

        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings. A larger value will lead to
          more information being printed.

        OUTPUT:

        A pAdicTree object that contatins that part of the given
        pAdicTree which satisfies the polynomial of this condition
        being equal to some power, at least least_exponent, of some
        number.

        If complement is set to True will also give the given tree as
        a second return value.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        T1, T0 = PolynomialCondition.pAdic_tree(self, pAdic_tree=pAdic_tree,
                                                pAdics=pAdics, complement=True,
                                                verbose=verbose, precision=1,
                                                **kwds)
        Te = PolynomialCondition.pAdic_tree(self, pAdic_tree=pAdic_tree,
                                            pAdics=pAdics, complement=False,
                                            verbose=verbose,
                                            precision=self.least_exponent(),
                                            **kwds)
        T = T0.union(Te)
        if complement:
            return T, pAdic_tree
        else:
            return T
    
    @cached_method
    def _x_str(self):
        i = 0
        x = 'x' + str(i)
        while x in self.variables():
            i += 1
            x = 'x' + str(i)
        return x

    @cached_method
    def _n_str(self):
        i = 0
        n = 'n' + str(i)
        while n in self.variables():
            i += 1
            n = 'n' + str(i)
        return n
        
    def _repr_(self):
        return ("The condition that " + str(self.polynomial()) + " == " +
                str(self._x_str()) + "^" + str(self._n_str()) + " with " +
                str(self._n_str()) + " >= " + str(self.least_exponent()))

    def _repr_short(self):
        return "%s == %s^%s"%(self.polynomial(), self._x_str(), self._n_str())
        
    def _latex_(self):
        return (latex(self.polynomial()) + " = " +
                "x_{" + self._x_str()[1:] + "}" +
                "^{n_{" + self._n_str()[1:] + "}}")

    def _cache_key(self):
        return 'PowerCondition', self.polynomial(), self.least_exponent()

    def __eq__(self, other):
        return (PolynomialCondition.__eq__(self, other) and
                self.least_exponent() == other.least_exponent())

class OrderCondition(PolynomialCondition):
    r"""The condition that the order of a polynomial is at most a given
    value at each prime.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import OrderCondition
        sage: R.<x, y> = ZZ[]
        sage: OrderCondition(x*y + y, n=2)
        The condition that x*y + y has order at most 2 at each prime.

    """

    def __init__(self, polynomial, n=0):
        r"""Initialize an OrderCondition.

        INPUT:

        - ``polynomial`` -- A polynomial for which the possible
          values should have order at most `n` at each prime.
        
        - ``n`` -- An integer (default=0) giving the maximal order
          `polynomial` can have at each prime.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import OrderCondition
            sage: R.<x, y> = ZZ[]
            sage: OrderCondition(x*y + y, n=2)
            The condition that x*y + y has order at most 2 at each prime.

        """
        PolynomialCondition.__init__(self, polynomial)
        self._n = n

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   verbose=False, **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those values for
        the variables such that the polynomial of this condition has
        order at most `n` at the prime of the corresponding p-adics.

        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings. A larger value will lead to
          more information being printed.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree which satisfies the polynomial of this condition
        having order at most `n` at the prime of the corresponding
        p-adics.

        If complement is set to True will also give the complement of
        the given Tree as a second return value.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        TN, TY = PolynomialCondition.pAdic_tree(self,
                                                pAdic_tree=pAdic_tree,
                                                pAdics=pAdics,
                                                complement=True,
                                                verbose=verbose,
                                                precision=self._n+1,
                                                **kwds)
        if complement:
            return TY, TN
        else:
            return TY
        
    def _repr_(self):
        return ("The condition that " + str(self.polynomial()) +
                " has order at most " + str(self._n) +
                " at each prime.") 

    def _repr_short(self):
        return "order of %s at most %s"%(self.polynomial(), self._n)
        
    def _latex_(self):
        return ("\\text{ord} \\left( " +
                latex(self.polynomial()) +
                "\\right) \\le " + latex(self._n))

    def _cache_key(self):
        return 'OrderCondition', self.polynomial(), self._n

    def __eq__(self, other):
        return (isinstance(other, OrderCondition) and
                PolynomialCondition.__eq__(self, other) and
                other._n == self._n)
    
class SquarefreeCondition(PolynomialCondition):
    r"""The condition that a polynomial is square free.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import SquarefreeCondition
        sage: R.<x, y> = ZZ[]
        sage: SquarefreeCondition(x*y + y)
        The condition that x*y + y is square free

    """

    def __init__(self, polynomial):
        r"""Initialize a SquarefreeCondition.

        INPUT:

        - ``polynomial`` -- A polynomial that should be an unknown
          square free number.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import SquarefreeCondition
            sage: R.<x, y> = ZZ[]
            sage: SquarefreeCondition(x*y + y)
            The condition that x*y + y is square free

        """
        PolynomialCondition.__init__(self, polynomial)

    def order(self):
        r"""Give the order of this Condition

        OUTPUT:

        An integer `n` for which the polynomial of this Condition has
        order at most `n` at each prime.

        """
        return self._n

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   verbose=False, **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those values for
        the variables such that the polynomial of this condition could
        be square free at those values.

        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings. A larger value will lead to
          more information being printed.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree which satisfies the polynomial of this condition
        being square free.

        If complement is set to True will also give the given tree as
        a second return value.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        TN, TY = PolynomialCondition.pAdic_tree(self, pAdic_tree=pAdic_tree,
                                                pAdics=pAdics, complement=True,
                                                verbose=verbose, precision=2,
                                                **kwds)
        if complement:
            return TY, pAdic_tree
        else:
            return TY
        
    def _repr_(self):
        return ("The condition that " + str(self.polynomial()) +
                " is square free")

    def _repr_short(self):
        return "%s square free"%(self.polynomial())
        
    def _latex_(self):
        return (latex(self.polynomial()) +
                " \\text{ is square free}")

    def _cache_key(self):
        return 'SquarefreeCondition', self.polynomial()

    def __eq__(self, other):
        return (isinstance(other, SquarefreeCondition) and
                PolynomialCondition.__eq__(self, other))

class CoprimeCondition(Condition_base):
    r"""The condition that variables are n-wise coprime.

    EXAMPLES::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: R.<x, y> = ZZ[]
        sage: CoprimeCondition([x,y])
        The condition that the variables ('x', 'y') are pairwise coprime.

    By default variables are assumed to be pairwise coprime, but other
    options are possible::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: R.<a, b, c, d> = ZZ[]
        sage: CoprimeCondition([a, b, c, d], n=0)
        The condition that always holds
        sage: CoprimeCondition([a, b, c, d], n=1)
        The condition that the variables ('a', 'b', 'c', 'd') are units.
        sage: CoprimeCondition([a, b, c, d], n=3)
        The condition that the variables ('a', 'b', 'c', 'd') are 3-wise coprime.
    """

    def __init__(self, variables, n=2):
        r""" The constructor of a CoprimeCondition.

        INPUT:
        
        - ``variables`` -- A collection of variables on which this
          condition applies. This may be any form of a variable, but
          will be converted into strings. Multiple variables with the
          same name are therefore not very well supported and may
          cause unpredictable behavior.

        - ``n`` -- A non-negative integer (default: 2) indicating the
          size of subsets of the variables that should be coprime,
          e.g. n=2 means that the variables should be pairwise coprime
          and n=1 indicates all variables should be units.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<x, y> = ZZ[]
            sage: CoprimeCondition([x,y])
            The condition that the variables ('x', 'y') are pairwise coprime.

        By default variables are assumed to be pairwise coprime, but
        other options are possible::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b, c, d> = ZZ[]
            sage: CoprimeCondition([a, b, c, d], n=0)
            The condition that always holds
            sage: CoprimeCondition([a, b, c, d], n=1)
            The condition that the variables ('a', 'b', 'c', 'd') are units.
            sage: CoprimeCondition([a, b, c, d], n=3)
            The condition that the variables ('a', 'b', 'c', 'd') are 3-wise coprime.
        """
        Condition_base.__init__(self, variables)
        self._n = n

    def number_of_coprimes(self):
        r"""Give the n for which the variables are n-wise coprime.

        OUTPUT:

        An integer n such that this condition is that the variables
        should be n-wise coprime.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b, c, d> = ZZ[]
            sage: CoprimeCondition([a, b, c, d], n=0).number_of_coprimes()
            0
            sage: CoprimeCondition([a, b, c, d], n=1).number_of_coprimes()
            1
            sage: CoprimeCondition([a, b, c, d], n=3).number_of_coprimes()
            3
            sage: CoprimeCondition([a, b, c, d]).number_of_coprimes()
            2
        """
        return self._n

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r""" Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree such that all
        variables in this condition are n-wise coprime with respect to
        the given prime.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree such that the variables in this condition are n-wise
        coprime.

        If complement was set to True will return a tuple with the
        afore mentioned as its first entry and the given pAdicTree as
        its second argument.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics,
                                   full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        T = pAdic_tree.root()
        tree_vars = pAdic_tree.variables()
        indices = tuple(tree_vars.index(var)
                        for var in self.variables()
                        if var in tree_vars)
        for node in T.children_at_level(1):
            if (sum(c == 0 for i,c in enumerate(node.quotient_tuple())
                    if i in indices) >= self._n):
                node.remove()
        if complement:
            return (pAdicTree(variables=pAdic_tree.variables(), root=T),
                    pAdic_tree)
        else:
            return pAdicTree(variables=pAdic_tree.variables(), root=T)

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        return self._n == 0

    def _repr_(self):
        if self._n == 0:
            return "The condition that always holds"
        if self._n == 1:
            return ("The condition that the variables " +
                    str(self.variables()) + " are units.")
        if self._n == 2:
            return ("The condition that the variables " +
                    str(self.variables()) + " are pairwise coprime.")
        return ("The condition that the variables " + str(self.variables()) +
                " are " + str(self._n) + "-wise coprime.")

    def _repr_short(self):
        if self._n == 0:
            return "true"
        if self._n == 1:
            return "%s are units"%(self.variables(),)
        if self._n == 2:
            return "%s are pairwise coprime"%(self.variables(),)
        return "%s are %s-wise coprime"%(self.variables(), self._n)
        
    def _latex_(self):
        if self._n == 0:
            return "\\top"
        if self._n == 1:
            return latex(self.variables()) + "\\text{ are units}"
        if self._n == 2:
            return latex(self.variables()) + "\\text{ are pairwise coprime}"
        return (latex(self.variables()) + "\\text{ are $" + str(self._n) +
                "$-wise coprime.}")

    def _cache_key(self):
        return 'CoprimeCondition', self.variables(), self._n

    def __eq__(self, other):
        return (Condition_base.__eq__(self, other) and
                self.number_of_coprimes() == other.number_of_coprimes())

class NotCondition(Condition_base):
    r""" The condition that another condition does not hold.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: R.<x, y> = ZZ[]
        sage: C = CoprimeCondition([x, y]); C
        The condition that the variables ('x', 'y') are pairwise coprime.
        sage: ~C
        The condition that the variables ('x', 'y') are not pairwise coprime.
    """

    def __init__(self, other):
        r"""Initializes a NotCondition.

        INPUT:

        - ``other`` -- The condition that this condition should be the
          negation of.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]); C
            The condition that the variables ('x', 'y') are pairwise coprime.
            sage: ~C
            The condition that the variables ('x', 'y') are not pairwise coprime.

        """
        self._other = other
        Condition_base.__init__(self, other.variables())

    def negated_condition(self):
        r"""Give the condition of which self is the negation

        OUTPUT:

        The condition that this condition is the negation of.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]); C
            The condition that the variables ('x', 'y') are pairwise coprime.
            sage: (~C).negated_condition()
            The condition that the variables ('x', 'y') are pairwise coprime.
        """
        return self._other

    def pAdic_tree(self, complement=False, **kwds):
        r""" Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of values for which
        the condition that this condition negates might not hold.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree which can possibly not satisfy the condition that
        this condition negates.

        If complement was set to True will return a tuple with the
        afore mentioned as its first entry and the tree of values for
        which the negated condition might hold.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        TY, TN = self._other.pAdic_tree(complement=True, **kwds)
        if complement:
            return TN, TY
        else:
            return TN

    def never(self):
        r"""Tell if this condition never holds

        .. NOTE::

        This function returning False does not imply this condition
        may not hold, as complex conditions might not be able to
        determine whether they never hold or not.

        OUTPUT:

        True or False. If the return value is True this condition can
        not hold on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.never()
            True
            sage: C2.never()
            False

        Note that when combining conditions, a condition that never
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that never holds
            sage: C1 | C2
            The condition that the variables ('a', 'b') are pairwise coprime.

        This method is not able to know of every combination of
        conditions whether it never holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 & C2).never()
            False

        """
        return self._other.always()

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        return self._other.never()

    def _repr_(self):
        s = self._other._repr_()
        s = s.replace(' are ', ' are not ')
        s = s.replace(' is ', ' is not ')
        s = s.replace(' not not ', ' ')
        s = s.replace(' always ', '<tmp>')
        s = s.replace(' never ', ' always ')
        s = s.replace('<tmp>', ' never ')
        s = s.replace('==', '<tmp>')
        s = s.replace('~=', '==')
        s = s.replace('<tmp>', '~=')
        s = s.replace(' and ', '<tmp>')
        s = s.replace(' or ', ' and ')
        s = s.replace('<tmp>', ' or ')
        return s

    def _repr_short(self):
        s = self._other._repr_short()
        s = s.replace(' are ', ' are not ')
        s = s.replace(' is ', ' is not ')
        s = s.replace(' not not ', ' ')
        s = s.replace(' always ', '<tmp>')
        s = s.replace(' never ', ' always ')
        s = s.replace('<tmp>', ' never ')
        s = s.replace('true', '<tmp>')
        s = s.replace('false', 'true')
        s = s.replace('<tmp>', 'false')
        s = s.replace('==', '<tmp>')
        s = s.replace('~=', '==')
        s = s.replace('<tmp>', '~=')
        s = s.replace(' and ', '<tmp>')
        s = s.replace(' or ', ' and ')
        s = s.replace('<tmp>', ' or ')
        return s
        
    def _latex_(self):
        s = self._other._latex_()
        s = s.replace(' are ', ' are not ')
        s = s.replace(' is ', ' is not ')
        s = s.replace(' not not ', ' ')
        s = s.replace('=', '<tmp>')
        s = s.replace('\\neq', '=')
        s = s.replace('<tmp>', '\\neq')
        s = s.replace('\\wedge', '<tmp>')
        s = s.replace('\\vee', '\\wedge')
        s = s.replace('<tmp>', '\\vee')
        s = s.replace('\\top', '<tmp>')
        s = s.replace('\\bot', '\\top')
        s = s.replace('<tmp>', '\\bot')
        return s

    def _cache_key(self):
        return 'NotCondition', self._other

    def __eq__(self, other):
        return (Condition_base.__eq__(self, other) and
                self.negated_condition() == other.negated_condition())

class AndCondition(Condition_base):
    r"""The condition that two conditions both hold.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
        sage: R.<x, y> = ZZ[]
        sage: CoprimeCondition([x,y]) & PolynomialCondition(x^2 + y^2 - 4)
        The condition that the variables ('x', 'y') are pairwise coprime and the condition that x^2 + y^2 - 4 == 0
    """

    def __init__(self, left, right):
        r""" Initialize an AndCondition.

        INPUT:

        - ``left`` -- The first condition that should hold for this
          condition to hold.

        - ``right`` - The second condition that should hold for this
          condition to hold.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: CoprimeCondition([x,y]) & PolynomialCondition(x^2 + y^2 - 4)
            The condition that the variables ('x', 'y') are pairwise coprime and the condition that x^2 + y^2 - 4 == 0
        """
        self._left = left
        self._right = right
        variables = list(left.variables())
        for var in right.variables():
            if var not in variables:
                variables.append(var)
        Condition_base.__init__(self, variables)

    def other_conditions(self):
        r"""Give the conditions that should hold for this condition to hold.

        OUTPUT:

        The conditions that should hold for this condition to hold.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C1 = CoprimeCondition([x, y]); C1
            The condition that the variables ('x', 'y') are pairwise coprime.
            sage: C2 = PolynomialCondition(x^2 + y^2 - 4); C2
            The condition that x^2 + y^2 - 4 == 0
            sage: (C1 & C2).other_conditions()
            (The condition that the variables ('x', 'y') are pairwise coprime.,
             The condition that x^2 + y^2 - 4 == 0)
        """
        return self._left, self._right

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a pAdicTree, returns the subtree that can satisfy both
        conditions this condition combines.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree which can satisfies both conditions combined in this
        condition.

        If complement was set to True will return a tuple with the
        afore mentioned as its first entry and as a second entry a
        pAdicTre object containing that part of the given pAdicTree
        for which either condition combined in this condition could
        not hold.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True
        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        T1 = self._left.pAdic_tree(pAdic_tree=pAdic_tree, pAdics=pAdics,
                                   complement=complement, **kwds)
        T2 = self._right.pAdic_tree(pAdic_tree=pAdic_tree, pAdics=pAdics,
                                    complement=complement, **kwds)
        if complement:
            return T1[0].intersection(T2[0]), T1[1].union(T2[1])
        else:
            return T1.intersection(T2)

    def never(self):
        r"""Tell if this condition never holds

        .. NOTE::

        This function returning False does not imply this condition
        may not hold, as complex conditions might not be able to
        determine whether they never hold or not.

        OUTPUT:

        True or False. If the return value is True this condition can
        not hold on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.never()
            True
            sage: C2.never()
            False

        Note that when combining conditions, a condition that never
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that never holds
            sage: C1 | C2
            The condition that the variables ('a', 'b') are pairwise coprime.

        This method is not able to know of every combination of
        conditions whether it never holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 & C2).never()
            False

        """
        return (self._left.never() or self._right.never())

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        return (self._left.always() and self._right.always())
        
    def _repr_(self):
        right_str = self._right._repr_()
        right_str = right_str[0].lower() + right_str[1:]
        left_str = self._left._repr_()
        if left_str.endswith('.'):
            left_str = left_str[:-1]
        return left_str + " and " + right_str

    def _repr_short(self):
        return self._left._repr_short() + " and " + self._right._repr_short()
        
    def _latex_(self):
        return self._left._latex_() + " \\wedge " + self._right._latex_()

    def _cache_key(self):
        return 'AndCondition', self._left, self._right

    def __eq__(self, other):
        return (Condition_base.__eq__(self, other) and
                self.other_conditions() == other.other_conditions())

class OrCondition(Condition_base):
    r"""The condition that either one of two conditions holds.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
        sage: R.<x, y> = ZZ[]
        sage: CoprimeCondition([x,y]) | PolynomialCondition(x^2 + y^2 - 4)
        The condition that the variables ('x', 'y') are pairwise coprime or the condition that x^2 + y^2 - 4 == 0
    """

    def __init__(self, left, right):
        r"""Initialize an OrCondition.

        INPUT:

        - ``left`` -- The first condition that could hold for this
          condition to hold.

        - ``right`` - The second condition that could hold for this
          condition to hold.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: CoprimeCondition([x,y]) | PolynomialCondition(x^2 + y^2 - 4)
            The condition that the variables ('x', 'y') are pairwise coprime or the condition that x^2 + y^2 - 4 == 0
        """
        self._left = left
        self._right = right
        variables = list(left.variables())
        for var in right.variables():
            if var not in variables:
                variables.append(var)
        Condition_base.__init__(self, variables)

    def other_conditions(self):
        r"""Give the conditions that can hold for this condition to hold.

        OUTPUT:

        The conditions that could hold for this condition to hold.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C1 = CoprimeCondition([x, y]); C1
            The condition that the variables ('x', 'y') are pairwise coprime.
            sage: C2 = PolynomialCondition(x^2 + y^2 - 4); C2
            The condition that x^2 + y^2 - 4 == 0
            sage: (C1 | C2).other_conditions()
            (The condition that the variables ('x', 'y') are pairwise coprime.,
             The condition that x^2 + y^2 - 4 == 0)
        """
        return self._left, self._right

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""Give this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of all those values
        for which at least one of the conditions combined in this
        condition can hold.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the given
        pAdicTree of which the values can satisfy at least one of the
        two conditions defined in this condition.

        If complement was set to True will return a tuple with the
        afore mentioned as its first entry and as the second entry a
        pAdicTree object that contains all the values of the given
        pAdicTree for which both condition might not hold.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        .. SEEALSO::

            :class:`pAdicTree`

        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        T1 = self._left.pAdic_tree(pAdic_tree=pAdic_tree, pAdics=pAdics,
                                   complement=complement, **kwds)
        T2 = self._right.pAdic_tree(pAdic_tree=pAdic_tree, pAdics=pAdics,
                                    complement=complement, **kwds)
        if complement:
            return T1[0].union(T2[0]), T1[1].intersection(T2[1])
        else:
            return T1.union(T2)

    def never(self):
        r"""Tell if this condition never holds

        .. NOTE::

        This function returning False does not imply this condition
        may not hold, as complex conditions might not be able to
        determine whether they never hold or not.

        OUTPUT:

        True or False. If the return value is True this condition can
        not hold on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.never()
            True
            sage: C2.never()
            False

        Note that when combining conditions, a condition that never
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that never holds
            sage: C1 | C2
            The condition that the variables ('a', 'b') are pairwise coprime.

        This method is not able to know of every combination of
        conditions whether it never holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 & C2).never()
            False

        """
        return (self._left.never() and self._right.never())

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        return (self._left.always() or self._right.always())
        
    def _repr_(self):
        right_str = self._right._repr_()
        right_str = right_str[0].lower() + right_str[1:]
        left_str = self._left._repr_()
        if left_str.endswith('.'):
            left_str = left_str[:-1]
        return left_str + " or " + right_str

    def _repr_short(self):
        return self._left._repr_short() + " or " + self._right._repr_short()
        
    def _latex_(self):
        return self._left._latex_() + " \\vee " + self._right._latex_()

    def _cache_key(self):
        return 'OrCondition', self._left, self._right

    def __eq__(self, other):
        return (Condition_base.__eq__(self, other) and
                self.other_conditions() == other.other_conditions())

class TreeCondition(Condition_base):
    r"""A condition that the values should be part of some pAdicTree.

    EXAMPLE::

        sage: from modular_method.padics.pAdic_base import pAdicBase
        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
        sage: R.<x, y> = ZZ[]
        sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 3))
        sage: TreeCondition(T)
        The condition that ('x', 'y') == (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2) mod 3

    """

    def __init__(self, pAdic_tree):
        r"""Initialize a TreeCondition.

        INPUT:

        - ``pAdic_tree`` -- A pAdicTree that contains the variables
          and values they should attain for this condition to hold.

        EXAMPLE::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 3))
            sage: TreeCondition(T)
            The condition that ('x', 'y') == (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2) mod 3

        """
        self._T = pAdic_tree
        Condition_base.__init__(self, pAdic_tree.variables())

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""Give this condition as a pAdicTree.

        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead. If that is also set
          to None, will use the pAdics of the tree stored in this
          Condition instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        OUTPUT:

        If no pAdicTree was given and no pAdics were given, returns
        the pAdicTree that defines this Condition. If complement was
        set to True will return that pAdicTree and its complement.

        If the given pAdicTree has no common pAdics with the pAdicTree
        stored in this Condition will return the given pAdicTree. If
        complement was set to True will return that pAdicTree twice.

        If the given pAdicTree has common pAdics with the pAdicTree
        stored in this Condition will return a pAdicTree containing
        all values of the given pAdicTree that agree with a value in
        the pAdicTree that defines this condition. Here two values of
        two pAdicTrees agree if the variables with the same name are
        assigned the same value.

        In the last case, if complement is set to True, will given the
        afore mentioned as the first return value and will give as the
        second return value a pAdicTree containing all values of the
        given pAdicTree that agree with a value of the complement of
        the pAdicTree that defines this condition.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        """
        if pAdic_tree is None:
            if pAdics is None:
                if complement:
                    return self._T, self._T.complement()
                else:
                    return self._T
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics,
                                   full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        if pAdics == self._T.pAdics():
            T1 = self._T.intersection(pAdic_tree)
            if complement:
                T2 = self._T.complement().intersection(pAdic_tree)
                return (T1.change_variables_to(pAdic_tree.variables()),
                        T2.change_variables_to(pAdic_tree.variables()))
            else:
                return T1.change_variables_to(pAdic_tree.variables())
        elif complement:
            return pAdic_tree, pAdic_tree
        else:
            return pAdic_tree

    def never(self):
        r"""Tell if this condition never holds

        .. NOTE::

        This function returning False does not imply this condition
        may not hold, as complex conditions might not be able to
        determine whether they never hold or not.

        OUTPUT:

        True or False. If the return value is True this condition can
        not hold on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.never()
            True
            sage: C2.never()
            False

        Note that when combining conditions, a condition that never
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = ~CoprimeCondition([a, b], n=0); C1
            The condition that never holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that never holds
            sage: C1 | C2
            The condition that the variables ('a', 'b') are pairwise coprime.

        This method is not able to know of every combination of
        conditions whether it never holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 & C2).never()
            False

        """
        return self._T.is_empty()

    def always(self):
        r"""Tell if this condition always holds

        .. NOTE::

        This function returning False does not imply this condition
        may have cases in which it does not hold, as complex
        conditions might not be able to determine whether they never
        hold or not.

        OUTPUT:

        True or False. If the return value is True this condition
        holds on any value for the variables.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1.always()
            True
            sage: C2.always()
            False

        Note that when combining conditions, a condition that always
        holds might make resulting expressions simpler::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: R.<a, b> = QQ[]
            sage: C1 = CoprimeCondition([a, b], n=0); C1
            The condition that always holds
            sage: C2 = CoprimeCondition([a, b]); C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 & C2
            The condition that the variables ('a', 'b') are pairwise coprime.
            sage: C1 | C2
            The condition that always holds

        This method is not able to know of every combination of
        conditions whether it always holds::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = QQ[]
            sage: C1 = CongruenceCondition(x, 2); C1
            The condition that x == 0 modulo 2
            sage: C2 = CongruenceCondition(x - 1, 2); C2
            The condition that x - 1 == 0 modulo 2
            sage: (C1 | C2).always()
            False

        """
        return self._T.is_full()

    def _repr_len(self, max_item=50, max_char=1000):
        r"""Give a string representation of this Condition of at most a given
        length

        INPUT:

        - ``max_item`` -- A non-negative integer (default: 50)
          indicating the maximal number of items to be included in
          this representation.

        - ``max_char`` -- A non-negative integer (default: 200) giving
          the maximal number of characters to be used in the string
          representation of this object.

        OUTPUT:

        The string representation of this condition or a shorter
        version thereof if that representation would have more than
        `max_char` characters or `max_item` items.

        """
        l = len(self.variables())
        if l == 0:
            return "true"
        values, modulus = self._T.give_as_congruence_condition()
        if hasattr(modulus, 'is_principal') and modulus.is_principal():
            modulus = modulus.gens_reduced()[0]
        if len(values) <= max_item:
            result = (str(self.variables() if l > 1 else self.variables()[0]) +
                      " == ")
            for i, value in enumerate(values):
                if len(result) > max_char:
                    break
                if i > 0:
                    result += ", "
                result += str(value if l > 1 else value[0])
            result += (" mod " +
                       (modulus._repr_short() if hasattr(modulus, '_repr_short')
                        else str(modulus)))
            if len(result) <= max_char:
                return result
        return (str(self.variables() if l > 1 else self.variables()[0]) +
                " is 1 of " + str(len(values)) + " possibilities mod " +
                (modulus._repr_short() if hasattr(modulus, '_repr_short')
                 else str(modulus)))
        
    def _repr_(self):
        result = "The condition that " + self._repr_len()
        result.replace("true", "always holds")
        result.replace("mod", "modulo")
        return result
               
    def _repr_short(self):
        return self._repr_len(max_item=10, max_char=40)

    def _latex_(self):
        l = len(self.variables())
        if l == 0:
            return "\\top"
        values, modulus = self._T.give_as_congruence_condition
        if hasattr(modulus, 'is_principal') and modulus.is_principal():
            modulus = modulus.gens_reduced()[0]
        result = (latex(self.variables() if l > 1 else self.variables()[0]) +
                  " = ")
        for i, value in enumerate(values):
            if i > 0:
                result += ", "
            result += latex(value if l > 1 else value[0])
        result += (" \\text{ (mod }" + latex(modulus) + "\\text{)}")
        return result

    def _cache_key(self):
        return 'TreeCondition', self._T

    def __and__(self, other):
        r"""Create the condition that both conditions hold.

        INPUT:
        
        - ``other`` -- A Condition, i.e. an instance of
          Condition_base.

        OUTPUT:

        A Condition object that holds on all values where both this
        Condition object and the given Condition object hold.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: CoprimeCondition([x,y]) & PolynomialCondition(x^2 + y^2 - 4)
            The condition that the variables ('x', 'y') are pairwise coprime and the condition that x^2 + y^2 - 4 == 0

        If two conditions are the same, the 'and' of both of them is
        just the first::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = ZZ[]
            sage: C1 = CongruenceCondition(x, 3); C1
            The condition that x == 0 modulo 3
            sage: C2 = CongruenceCondition(x, 3); C2
            The condition that x == 0 modulo 3
            sage: C1 & C2
            The condition that x == 0 modulo 3

        .. SEEALSO ::
        
           :class:`AndCondition`

        """
        if (isinstance(other, TreeCondition) and self._T.pAdics() == other._T.pAdics()):
            return TreeCondition(self._T.intersection(other._T))
        else:
            return Condition_base.__and__(self, other)

    def __or__(self, other):
        r"""Create the condition that either condition holds.

        INPUT:
        
        - ``other`` -- A Condition, i.e. an instance of
          Condition_base.

        OUTPUT:

        A Condition object that holds on all values where either this
        Condition object or the given Condition object holds.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: CoprimeCondition([x,y]) | PolynomialCondition(x^2 + y^2 - 4)
            The condition that the variables ('x', 'y') are pairwise coprime or the condition that x^2 + y^2 - 4 == 0

        If two conditions are the same, the 'or' of both of them is
        just the first::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x> = ZZ[]
            sage: C1 = CongruenceCondition(x, 3); C1
            The condition that x == 0 modulo 3
            sage: C2 = CongruenceCondition(x, 3); C2
            The condition that x == 0 modulo 3
            sage: C1 | C2
            The condition that x == 0 modulo 3

        .. SEEALSO::

            :class:`OrCondition`

        """
        if (isinstance(other, TreeCondition) and self._T.pAdics() == other._T.pAdics()):
            return TreeCondition(self._T.union(other._T))
        else:
            return Condition_base.__or__(self, other)

    def __invert__(self):
        r"""Create the condition that this condition does not hold.

        OUTPUT:

        A Condition object that holds on all values where this
        Condition does not hold.

        EXAMPLES::
        
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: ~ PolynomialCondition(x^2 + y^2 - 4)
            The condition that x^2 + y^2 - 4 ~= 0

        Note that a double not simplifies in print, but gives a
        different object::

            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1); C
            The condition that -x^3 + y^2 - 1 == 0
            sage: ~~C
            The condition that -x^3 + y^2 - 1 == 0
            sage: C == ~~C
            False

        .. SEEALSO::

            :class:`NotCondition`

        """
        return TreeCondition(self._T.complement())

    def __eq__(self, other):
        return (Condition_base.__eq__(self, other) and
                self.pAdic_tree() == other.pAdic_tree())

class TextCondition(Condition_base):
    r"""A condition described by text.

    EXAMPLES::

        sage: from modular_method.diophantine_equations.conditions import TextCondition
        sage: C = TextCondition("The condition that the ball is blue"); C
        The condition that the ball is blue
        sage: ~C
        The condition that the ball is not blue
        sage: C = TextCondition("The condition that n is a square", "n == 'square'"); C
        The condition that n is a square
        sage: C._repr_short()
        "n == 'square'"

    """

    def __init__(self, text, short_text=None, variables=[]):
        r"""Initialize a TextCondition

        INPUT:

        - ``text`` -- A string that is the string representation of
          this condition.
        
        - ``short_text`` -- A string or `None` (default) that is the
          short string representation of this condition. If set to
          `None` it will be set to be the same as `text`.
        
        - ``variables`` -- A list or tuple of variables (default: [])
          on which this condition applies. This is not necessary for
          this condition, but is implemented for compatibility.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import TextCondition
            sage: C = TextCondition("The condition that the ball is blue"); C
            The condition that the ball is blue
            sage: ~C
            The condition that the ball is not blue
            sage: C = TextCondition("The condition that n is a square", "n == 'square'"); C
            The condition that n is a square
            sage: C._repr_short()
            "n == 'square'"

        """
        self._long = text
        if short_text == None:
            self._short = self._long
        else:
            self._short = short_text
        Condition_base.__init__(self, variables)

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""Give this condition as a pAdicTree.

        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None) on which
          this condition should be applied. If set to None will be
          initiated as the full tree with the given pAdics.

        - ``pAdics`` -- A pAdicBase object (default: None) determining
          the pAdics that should be used. If set to None will use the
          pAdics of the given pAdicTree instead. If that is also set
          to None, will use the pAdics of the tree stored in this
          Condition instead.

        - ``complement`` -- A boolean (default: False) determining
          whether the complement of the result should be returned.

        OUTPUT:

        If no pAdicTree was given and no pAdics were given, returns
        the pAdicTree that defines this Condition. If complement was
        set to True will return that pAdicTree and its complement.

        If the given pAdicTree has no common pAdics with the pAdicTree
        stored in this Condition will return the given pAdicTree. If
        complement was set to True will return that pAdicTree twice.

        If the given pAdicTree has common pAdics with the pAdicTree
        stored in this Condition will return a pAdicTree containing
        all values of the given pAdicTree that agree with a value in
        the pAdicTree that defines this condition. Here two values of
        two pAdicTrees agree if the variables with the same name are
        assigned the same value.

        In the last case, if complement is set to True, will given the
        afore mentioned as the first return value and will give as the
        second return value a pAdicTree containing all values of the
        given pAdicTree that agree with a value of the complement of
        the pAdicTree that defines this condition.

        EXAMPLES::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y]) & PolynomialCondition(x^2 + y^2 - 4)
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 3), precision=3)
            sage: T.get_values_at_level(1)
            [(0, 1), (0, 2), (1, 0), (2, 0)]

        The complement can be used to get two sets, one for which the
        condition is satisfied and one for which it is not::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(y^2 - x^3 - 1)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True, precision=3)
            sage: Ty.get_values_at_level(1)
            [(0, 1), (1, 0)]
            sage: Tn.get_values_at_level(1)
            [(0, 0), (1, 0), (1, 1)]

        One can use custom trees to limit the values on which a
        condition should be applied::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition
            sage: R.<x, y> = ZZ[]
            sage: C = PolynomialCondition(x^2 + y^2 - 4)
            sage: C.pAdic_tree(pAdics=pAdicBase(QQ, 2), precision=2).get_values_at_level(1)
            [(0, 0)]
            sage: T = CoprimeCondition([x, y]).pAdic_tree(pAdics=pAdicBase(QQ, 2))
            sage: C.pAdic_tree(pAdic_tree=T, precision=2).get_values_at_level(1)
            []

        Some Condition objects accept that both the pAdic_tree
        argument and pAdics argument are set to None, but only in case
        it is obvious which tree should be returned::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, TreeCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: T = C.pAdic_tree(pAdics=pAdicBase(QQ, 5))
            sage: C2 = TreeCondition(T)
            sage: C2.pAdic_tree()
            p-adic tree for the variables ('x', 'y') with respect to p-adics given by Rational Field and (5)
            sage: C.pAdic_tree()
            Traceback (most recent call last):
            ...
            ValueError: At least the argument prime must be set

        The complement returned might not in all cases be disjoint
        from the first tree::

            sage: from modular_method.padics.pAdic_base import pAdicBase
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x^2 + 2*y^2, 3)
            sage: Ty, Tn = C.pAdic_tree(pAdics=pAdicBase(QQ, 2), complement=True)
            sage: Ty == Tn
            True

        """
        # Only implemented for compatibility in elimination functions
        # Will return the entire given tree as result and as complement
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics, full=True)
        if complement:
            return pAdic_tree, pAdic_tree
        else:
            return pAdic_tree


    def _repr_(self):
        return self._long

    def _repr_short(self):
        return self._short

    def __eq__(self, other):
        return (isinstance(other, TextCondition) and
                isinstance(self, other.__class__) and
                self._long == other._long and
                self._short == other._short)

class ConditionalValue(SageObject):
    r"""Some value that depends on some condition.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue
        sage: R.<x, y> = ZZ[]
        sage: C = CoprimeCondition([x, y])
        sage: ConditionalValue([(1, C), (0, ~C)])
        1 if ('x', 'y') are pairwise coprime
        0 if ('x', 'y') are not pairwise coprime

    """

    def __init__(self, val_con):
        r"""Initialize a ConditionalValue

        INPUT:

        - ``val_con`` - A list of tuples, where each tuple consists of
          a value and a condition on when this value is attained, in
          that order.  A value can be any object, whilst a condition
          must extend Condition_base. The different conditions do not
          have to include all possibilities, nor do they have to be
          exclusive of one another, but not adhering to this will
          reflect in the resulting object.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: ConditionalValue([(1, C), (0, ~C)])
            1 if ('x', 'y') are pairwise coprime
            0 if ('x', 'y') are not pairwise coprime

        """
        self._vals = tuple(vc[0] for vc in val_con)
        self._con = tuple(vc[1] for vc in val_con)
        for c in self._con:
            if not isinstance(c, Condition_base):
                raise ValueError("%s is not a condition"%(c,))

    def no_value(self):
        r"""Tell if this conditional value does not attain any value.

        OUTPUT:

        True - If this condition does not contain any possible value.

        False - If this condition contains at least one value.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x + y, 4)
            sage: V = ConditionalValue([(12, C)]); V
            12 if x + y == 0 mod 4
            sage: V.no_value()
            False
            sage: ConditionalValue([]).no_value()
            True

        """
        return len(self._vals) == 0

    def _repr_lines(self):
        if self.no_value():
            return []
        result = [str(val) for val in self._vals]
        l = max(len(r) for r in result) + 1
        result = [r + (' '*(l - len(r))) + "if " + self._con[i]._repr_short()
                  for i, r in enumerate(result)]
        return result
    
    def _repr_(self):
        if self.no_value():
            return "Conditional value with no possible value."
        lines = self._repr_lines()
        result = ""
        for i, line in enumerate(lines):
            if i > 0:
                result += "\n"
            result += line
        return result

    def _latex_lines(self):
        if self.no_value():
            return []
        return [latex(val) + "& \\text{ if }" + latex(self._con[i])
                for i, val in enumerate(self._vals)]
    
    def _latex_(self):
        if self.no_value():
            return "\\text{Conditional value with no possible value.}"
        result = "\\left\\{ \\begin{array}{lr}\n"
        for i, line in enumerate(self._latex_lines()):
            if i > 0:
                result += " \\\\\n"
            result += line
        result += "\n \\end{array} \\right."
        return result

    def __add__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR, self, other)

    def __radd__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR, other, self)

    def __sub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, self, other)

    def __rsub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, other, self)

    def __mul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, self, other)

    def __rmul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, other, self)

    def __div__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, self, other)

    def __rdiv__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, other, self)

    def __pow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, self, other)

    def __rpow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, other, self)

    def __iter__(self):
        return iter(zip(self._vals, self._con))

    def __len__(self):
        return len(self._vals)

    def __getitem__(self, index):
        return (self._vals[index], self._con[index])

    def _cache_key(self):
        return tuple((val, 'if', con) for val, con in self)

    def __eq__(self, other):
        return (isinstance(other, ConditionalValue) and
                len(self) == len(other) and
                all(self[i] == other[i] for i in range(len(self))))

    def __ne__(self, other):
        return not self.__eq__(other)

class ConditionalExpression(SageObject):
    SUM_OPERATOR = ('+', '+', 0)
    MINUS_OPERATOR = ('-', '-', 0.5)
    PRODUCT_OPERATOR = ('*', '\\cdot', 2)
    DIVISION_OPERATOR = ('/','/', 2.5)
    EXPONENT_OPERATOR = ('^', '^', 4.5)
    r"""An expression containing conditional values.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CongruenceCondition, ConditionalValue
        sage: R.<x, y> = ZZ[]
        sage: C = CongruenceCondition(x + y, 4)
        sage: V = ConditionalValue([(2, C), (-2, ~C)])
        sage: V * 3 + 12
        n0*3+12
         where 
        n0 = 2  if x + y == 0 mod 4
             -2 if x + y ~= 0 mod 4

    """
    
    def __init__(self, operator, left, right):
        r"""Initialize a conditional expression.

        INPUT:

        - ``operator`` -- A tuple containing in this order a string
          representing the operator, a string that will produce the
          operator in latex and a non-negative integer indicating the
          power of the operator.

        EXAMPLES::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x + y, 4)
            sage: V = ConditionalValue([(2, C), (-2, ~C)])
            sage: V * 3 + 12
            n0*3+12
             where 
            n0 = 2  if x + y == 0 mod 4
                 -2 if x + y ~= 0 mod 4

        ConditionalExpressions can be used to write formulas with
        objects that normally don't support this::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: ConditionalValue([(2, C),(0, ~C)]) + "I"
            n0+I
             where 
            n0 = 2 if ('x', 'y') are pairwise coprime
                 0 if ('x', 'y') are not pairwise coprime

        Constructing expressions explicitly allows for custom symbols
        to be used::

            sage: from modular_method.diophantine_equations.conditions import ConditionalExpression
            sage: E = ConditionalExpression((' & ', ' \\alpha ', 1), "Hello", "World"); E
            Hello & World
            sage: latex(E)
            Hello \alpha World

        """
        self._op = operator
        self._left = left
        self._right = right

    def left(self):
        r"""Give the left side of this expression.

        A ConditionalExpression consists of a left and a right side
        separated by some operator.

        OUTPUT:

        The left side of this expression.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: E = ConditionalValue([(2, C), (-2, ~C)]) + "I"; E
            n0+I
             where 
            n0 = 2  if ('x', 'y') are pairwise coprime
                 -2 if ('x', 'y') are not pairwise coprime
            sage: E.left()
            2  if ('x', 'y') are pairwise coprime
            -2 if ('x', 'y') are not pairwise coprime
            sage: E.right()
            'I'
            sage: E.operator()
            ('+', '+', 0)

        """
        return self._left

    def operator(self):
        r"""Give the operator of this expression.

        A ConditionalExpression consists of a left and a right side
        separated by some operator.

        OUTPUT:

        A tuple consisting of the string representing the operator in
        this expression in sage, the string representing the operator
        in LaTeX and the priority of this operator.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: E = ConditionalValue([(2, C), (-2, ~C)]) + "I"; E
            n0+I
             where 
            n0 = 2  if ('x', 'y') are pairwise coprime
                 -2 if ('x', 'y') are not pairwise coprime
            sage: E.left()
            2  if ('x', 'y') are pairwise coprime
            -2 if ('x', 'y') are not pairwise coprime
            sage: E.right()
            'I'
            sage: E.operator()
            ('+', '+', 0)

        """
        return self._op

    def right(self):
        r"""Give the right side of this expression.

        A ConditionalExpression consists of a left and a right side
        separated by some operator.

        OUTPUT:

        The right side of this expression.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CoprimeCondition([x, y])
            sage: E = ConditionalValue([(2, C), (-2, ~C)]) + "I"; E
            n0+I
             where 
            n0 = 2  if ('x', 'y') are pairwise coprime
                 -2 if ('x', 'y') are not pairwise coprime
            sage: E.left()
            2  if ('x', 'y') are pairwise coprime
            -2 if ('x', 'y') are not pairwise coprime
            sage: E.right()
            'I'
            sage: E.operator()
            ('+', '+', 0)

        """
        return self._right

    def _factor_side(self, side):
        r"""
        Factors a side. See :meth:`factors` for details.
        """
        if isinstance(side, ConditionalExpression):
            return side.factors()
        if hasattr(side, 'factor'):
            F = side.factor()
            result = {f: e for f,e in F}
            if hasattr(F, 'unit') and F.unit() != 1:
                result[F.unit()] = 1
            return result
        return {side: 1}
    
    def factors(self):
        r"""Give the factors in this expression and their exponents.

        The factoring will work as follows:

        - If the operator is a PRODUCT_OPERATOR will factor left and
          right and combine both by taking the union of the factors of
          left and right. The exponent of a factor is the sum of the
          exponents if the factor appeared on both sides or the
          exponent of the side it appeared on otherwise.

        - If the operator is a DIVISION_OPERATOR will factor left and
          right and combine both by taking the union of the factors of
          left and right. The exponent of a factor is exponent of that
          factor in left minus the exponent of that factor in
          right. If the factor appeared only on the left will take the
          corresponding exponent instead and if the factor only
          appeared on the right will take 0 minus the corresponding
          exponent instead.

        - If the operator is an EXPONENT_OPERATOR will factor left
          and take all those factors. The exponent of each factor
          will be the corresponding exponent in left times the
          right side.

        - If the operator is anything else, will return this as the
          only factor and 1 as the exponent.

        - When factoring a side, will use the method factors if the
          side is a ConditionalExpression, will attempt to use the
          method factor if it is any other object with such a method
          and in any other case will assume that side to be a factor
          with exponent 1. If the method factor is used and the
          resulting object has a method unit, the result of this
          method will be added as a factor with respective exponent 1
          unless it is 1 itself.

        OUTPUT:

        A dictionary containing as keys the factors of this expression
        and as values their corresponding exponents.  Note that these
        factors and exponents could be any type of expression that
        could normally appear in a conditional expression.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x + 2*y, 5)
            sage: V = ConditionalValue([(2, C), (3, ~C)])
            sage: E = 24^V
            sage: E.factors()
            {2: 3*n0
              where 
             n0 = 2 if x + 2*y == 0 mod 5
                  3 if x + 2*y ~= 0 mod 5, 3: 2 if x + 2*y == 0 mod 5
             3 if x + 2*y ~= 0 mod 5}

        """
        if self._op == ConditionalExpression.PRODUCT_OPERATOR:
            result = self._factor_side(self._left)
            extra = self._factor_side(self._right)
            for f in extra:
                if f in result and result[f] != 0:
                    result[f] = result[f] + extra[f]
                elif extra[f] != 0:
                    result[f] = extra[f]
            return result
        if self._op == ConditionalExpression.DIVISION_OPERATOR:
            result = self._factor_side(self._left)
            extra = self._factor_side(self._right)
            for f in extra:
                if f in result and result[f] != 0:
                    result[f] = result[f] - extra[f]
                elif extra[f] != 0:
                    result[f] = 0 - extra[f]
            return result
        if self._op == ConditionalExpression.EXPONENT_OPERATOR:
            result = self._factor_side(self._left)
            for f in result:
                if result[f] == 1:
                    result[f] = self._right
                elif self._right != 1:
                    result[f] = result[f] * self._right
            return result
        return {self: 1}

    def _repr_side(self, side, vals, bracket_level):
        if isinstance(side, ConditionalExpression):
            return side._repr_info(vals, bracket_level)
        if isinstance(side, ConditionalValue):
            vals.append(side)
            return "n" + str(len(vals) - 1)
        if hasattr(side, "_repr_short"):
            return side._repr_short()
        return str(side)

    def _operation_on(self, left, right):
        r"""Perform the operation of this expression on two values.

        Depending on the operation will do the following

        - If the operation is SUM_OPERATION will return the sum of
          left plus right

        - If the operation is MINUS_OPERATION will return left minus
          right.

        - If the operation is PRODUCT_OPERATION will return the
          product of left with right.

        - If the operation is DIVISION_OPERATION will return left
          divided by right.

        - If the operation is EXPONENT_OPERATION will return left to
          the power right

        - If the operation is anything else will return a ValueError.

        """
        if self._op == ConditionalExpression.SUM_OPERATOR:
            return left + right
        if self._op == ConditionalExpression.MINUS_OPERATOR:
            return left - right
        if self._op == ConditionalExpression.PRODUCT_OPERATOR:
            return left * right
        if self._op == ConditionalExpression.DIVISION_OPERATOR:
            return left / right
        if self._op == ConditionalExpression.EXPONENT_OPERATOR:
            return left ^ right
        raise ValueError('Can not evaluate operation %s'%(self._op,))
    
    def value(self):
        r"""Turn this ConditionalExpression into a ConditionalValue.

        Attempts to turn this ConditionalExpression into a
        ConditionalValue by taking for each ConditionalValue in the
        expression a list of possible combinations of their conditions
        as the conditions for the resulting conditional value. For
        each of these conditions will substitute the appropiate value
        for each ConditionalValue in this expression and evaluate the
        expression assigning the result as the corresponding value of
        the resulting ConditionalValue.

        .. NOTE::

        For this method to work all operations in this
        ConditionalExpression must be one of the predefined operations
        corresponding to addition, subtraction, multiplication,
        division and exponentiation, and the objects in the Expression
        must allow these operations to be performed on them.

        OUTPUT:

        A ConditionalValue or some element that has exactly the same
        values as this expression on every choice of parameters.

        EXAMPLE::

            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition, ConditionalValue
            sage: R.<x, y> = ZZ[]
            sage: C = CongruenceCondition(x + 2*y, 5)
            sage: V = ConditionalValue([(2, C), (3, ~C)])
            sage: E = 2 * V + 4; E
            2*n0+4
             where 
            n0 = 2 if x + 2*y == 0 mod 5
                 3 if x + 2*y ~= 0 mod 5
            sage: E.value()
            8  if x + 2*y == 0 mod 5
            10 if x + 2*y ~= 0 mod 5

        """
        left = self._left
        if isinstance(left, ConditionalExpression):
            left = left.value()
        right = self._right
        if isinstance(right, ConditionalExpression):
            right = right.value()
        if isinstance(left, ConditionalValue):
            result = []
            for val_l, con_l in left:
                if isinstance(right, ConditionalValue):
                    for val_r, con_r in right:
                        result.append((self._operation_on(val_l, val_r),
                                       con_l & con_r))
                else:
                    result.append((self._operation_on(val_l, right), con_l))
            return ConditionalValue(result)
        if isinstance(right, ConditionalValue):
            result = []
            for val, con in right:
                result.append((self._operation_on(left, val), con))
            return ConditionalValue(result)
        return self._operation_on(left, right)
        
    def _repr_info(self, vals, bracket_level):
        result = ""
        if bracket_level > self._op[2]:
            result += "("
        result += self._repr_side(self._left, vals, floor(self._op[2]))
        result += self._op[0]
        result += self._repr_side(self._right, vals, ceil(self._op[2]))
        if bracket_level > self._op[2]:
            result += ")"
        return result

    def _repr_(self):
        vals = []
        result = self._repr_info(vals, 0)
        if len(vals) > 0:
            result += "\n where \n"
            front_space = ceil(ZZ(len(vals)).log(10)) + 5
            for i, val in enumerate(vals):
                lines = val._repr_lines()
                if len(lines) == 0:
                    if i > 0:
                        result += "\n"
                    result += "n" + str(i) + " ="
                    result += " no possible value"
                for j, line in enumerate(lines):
                    if i + j > 0:
                        result += "\n"
                    r = ""
                    if j == 0:
                        r += "n" + str(i) + " ="
                    result += r + " "*(front_space-len(r)) + \
                              line
        return result

    def _latex_side(self, side, vals, bracket_level):
        if isinstance(side, ConditionalExpression):
            return side._latex_info(vals, bracket_level)
        if isinstance(side, ConditionalValue):
            vals.append(side)
            return "n_{" + str(len(vals) - 1) + "}"
        return str(side)
        
    def _latex_info(self, vals, bracket_level):
        result = ""
        if bracket_level > self._op[2]:
            result += "\\left("
        result += self._latex_side(self._left, vals, floor(self._op[2]))
        result += self._op[1]
        result += self._latex_side(self._right, vals, ceil(self._op[2]))
        if bracket_level > self._op[2]:
            result += "\\right)"
        return result

    def _latex_(self):
        vals = []
        result = self._latex_info(vals,0)
        if len(vals) > 0:
            result += "\n\\\\\\text{ where } \\\\\n"
            for i, val in enumerate(vals):
                if i > 0:
                    result += " \\\n"
                result += "n_{" + str(i+1) + "} = "
                result += " & " + latex(val)
        return result

    def __add__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR,self, other)

    def __radd__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR, other, self)

    def __sub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, self, other)

    def __rsub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, other, self)

    def __mul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, self, other)

    def __rmul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, other, self)

    def __div__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, self, other)

    def __rdiv__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, other, self)

    def __pow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, self, other)

    def __rpow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, other, self)

    def __eq__(self, other):
        return (isinstance(other, ConditionalExpression) and
                self.left() == other.left() and self.right() == other.right())

    def __ne__(self, other):
        return not self.__eq__(other)

def apply_to_conditional_value(function, value, singleton=False,
                               use_condition=False, default_condition=None):
    r"""Apply a function to a conditional value.

    INPUT:
    
    - ``function`` -- Any function with a single input or two inputs
      if the argument use_condition is set to True.

    - ``value`` -- Any value that is accepted by the given function as
      an input or a ConditionalValue that contains such values.

    - ``singleton`` -- A boolean value (default: False) indicating
      whether a ConditionalValue with only one possibility should be
      returned as a ConditionalValue. If set to False and the
      resulting ConditionalValue only has one value, will just return
      that value.

    - ``use_condition`` -- A boolean value (default: False). If set to
      true, will pass a value condition pair to the function for each
      possible value instead of only the value. This is useful for
      functions which also require the condition on which this value
      depends.

    - ``default_condition`` -- If the given value was not a
      ConditionalValue, will use this to determine the condition
      corresponding to the value if it is required to be passed to the
      funcion or to construct a ConditionalValue at the end.

    OUTPUT:

    The function evaluated on the given value.

    If the given value was a ConditionalValue this means the function
    is evaluated on every value in the ConditionalValue producing a
    new ConditionalValue of possible outcomes. Conditions which
    produce the same outcome will be combined using an OrCondition.

    If the resulting ConditionalValue would have only one possibility
    and singleton is set to False, will return the value of that
    single possibility instead of the whole ConditionalValue.

    If the return value of the given function at any of the tried
    values is itself a ConditionalValue `V`, the ConditionalValue
    returned by this function will contain each the values in `V` as a
    possible outcome rather than using `V` as an outcome. Note that
    the conditions for each of these outcomes will be the same as
    those in `V`.

    EXAMPLES::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue, apply_to_conditional_value
        sage: R.<x, y> = ZZ[]
        sage: C = CoprimeCondition([x, y])
        sage: V = ConditionalValue([(1, C), (0, ~C)])
        sage: def f(x):
        ....:     return 3 * x
        ....: 
        sage: apply_to_conditional_value(f, V)
        3 if ('x', 'y') are pairwise coprime
        0 if ('x', 'y') are not pairwise coprime

    Acts as the function on normal values::

        sage: from modular_method.diophantine_equations.conditions import apply_to_conditional_value
        sage: def f(x):
        ....:     return 3 * x
        ....: 
        sage: apply_to_conditional_value(f, 2)
        6

    The function can use the condition::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, ConditionalValue, apply_to_conditional_value
        sage: R.<x, y> = ZZ[]
        sage: C = CoprimeCondition([x, y])
        sage: V = ConditionalValue([(1, C), (0, ~C)])
        sage: def f(x, con):
        ....:     return (3*x if con == C else x - 8)
        ....: 
        sage: apply_to_conditional_value(f, V, use_condition=True)
        3  if ('x', 'y') are pairwise coprime
        -8 if ('x', 'y') are not pairwise coprime

    A single answer is by default just returned as the value, but can
    be forced to be a ConditionalValue::

        sage: from modular_method.diophantine_equations.conditions import CongruenceCondition, ConditionalValue, apply_to_conditional_value
        sage: R.<x, y> = ZZ[]
        sage: C = CongruenceCondition(x + 2*y, 5)
        sage: V = ConditionalValue([(-2, C), (2, ~C)])
        sage: def f(x):
        ....:     return x^2
        ....: 
        sage: apply_to_conditional_value(f, V)
        4
        sage: apply_to_conditional_value(f, V, singleton=True)
        4 if x + 2*y == 0 mod 5 or x + 2*y ~= 0 mod 5

    The default condition has to be given in some cases for the answer
    to be as expected::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, apply_to_conditional_value
        sage: R.<x, y> = ZZ[]
        sage: C = CoprimeCondition([x, y])
        sage: def f(x, con):
        ....:     return (6 - x if con == C else 0)
        ....: 
        sage: apply_to_conditional_value(f, 3, use_condition=True)
        0
        sage: apply_to_conditional_value(f, 3, use_condition=True, default_condition=C)
        3
        sage: def g(x):
        ....:     return 6 + x
        ....: 
        sage: apply_to_conditional_value(g, 3, singleton=True)
        Traceback (most recent call last):
        ...
        ValueError: None is not a condition
        sage: apply_to_conditional_value(g, 3, singleton=True, default_condition=C)
        9 if ('x', 'y') are pairwise coprime

    If the outcome of the provided function is a conditional value,
    the resulting conditional value will produce a single conditional
    value that incorporates that outcome::

        sage: from modular_method.diophantine_equations.conditions import CongruenceCondition, ConditionalValue, apply_to_conditional_value
        sage: R.<a, b> = ZZ[]
        sage: C1 = CongruenceCondition(a, 2)
        sage: C2 = CongruenceCondition(b, 2)
        sage: V1 = ConditionalValue([(1, C1), (2, ~C1)])
        sage: V2 = ConditionalValue([(0, C2), (1, ~C2)])
        sage: def f(x):
        ....:     return (V1 if x == 0 else 3)
        ....:
        sage: apply_to_conditional_value(f, V2)
        1 if a == 0 mod 2 
        2 if a ~= 0 mod 2
        3 if b ~= 0 mod 2

    .. SEEALSO::

        :class:`ConditionalValue`

    """
    if not isinstance(value, ConditionalValue):
        value = [(value, default_condition)]
    values = []
    conditions = []
    for val, con in value:
        if use_condition:
            result = function(val, con)
        else:
            result = function(val)
        if not isinstance(result, ConditionalValue):
            result = [(result, con)]
        for f_val, f_con in result:
            if (f_con is None or  not f_con.never()):
                try:
                    i = values.index(f_val)
                    conditions[i] = conditions[i] | f_con
                except ValueError:
                    values.append(f_val)
                    conditions.append(f_con)
    if not singleton and len(values) == 1:
        return values[0]
    else:
        return ConditionalValue(list(zip(values, conditions)))

def conditional_over_values(function, values, start_condition=None,
                            singleton=False):
    r"""Construct one conditional value by evaluating a function on given
    values.

    A function of which the outcome depends on some unknown parameters
    might have cases for the parameters in which it does not work. In
    such a case one might try to apply the function again to a
    different value. Iterating this process over different values will
    produce different return values for many different cases. This
    function will do this and return all these values in a single
    ConditionalValue with the corresponding cases.

    INPUT:

    - ``function`` -- A function with two inputs, a value and a
      condition, which has a single return value. The return value
      must be the possible return values of this function at a given
      value if the given condition is satisfied. If the function would
      give no return value, the return value `None` may be used. If
      multiple return values are possible, all of them should be
      combined in a ConditionalValue giving these values and
      conditions when they occur.

    - ``values`` -- An iterable object consisting of values. These may
      not be ConditionalValues.

    - ``start_condition`` -- A Condition or `None` (default) which
      should be the condition passed to the given function when
      evaluated on the first value of `values`. Note that this can
      only be `None` if the given function accepts `None` as a
      condition.

    - ``singleton`` -- A boolean value (default: `False`) which
      indicates whether a single value should still be returned as a
      ConditionalValue. If set to False a ConditionalValue with only a
      single value will just be returned as that single value.

    OUTPUT:

    A ConditionalValue consisting of values `V` with corresponding
    conditions `C` such that `V` is the first return value that is not
    `None` when `function` is applied to each value in `values` under
    the condition `C`. A value `V` may only be `None` if `function`
    returns `None` for each value in `values` under the condition
    `C`.

    EXAMPLE::

        sage: from modular_method.padics.pAdic_base import pAdicBase
        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
        sage: from modular_method.diophantine_equations.conditions import ConditionalValue
        sage: from modular_method.diophantine_equations.conditions import TreeCondition
        sage: from modular_method.diophantine_equations.conditions import conditional_over_values
        sage: R.<x> = ZZ[]
        sage: C = PolynomialCondition(x^2 + 1)
        sage: pAdics = pAdicBase(QQ, 5)
        sage: def f(val, con):
        ....:     Y, N = C.pAdic_tree(pAdic_tree=con.pAdic_tree(pAdics=pAdics),
        ....:                         complement=True, precision=val)
        ....:     success = "Solution for precision " + str(val)
        ....:     return ConditionalValue([(success, TreeCondition(Y)),
        ....:                              (None, TreeCondition(N))])
        ....: 
        sage: start = CoprimeCondition([x])
        sage: conditional_over_values(f, [4, 3, 2, 1], start_condition=start)
        Solution for precision 4 if x == 182, 443 mod 625
        Solution for precision 3 if x is 1 of 8 possibilities mod 625
        Solution for precision 2 if x is 1 of 8 possibilities mod 125
        Solution for precision 1 if x == 2, 3, 8, 12, 13, 17, 22, 23 mod 25
        None                     if x == 0, 1, 4 mod 5

    Setting the start_condition to `None` only works if the given
    function supports `None` as a condition argument::

        sage: from modular_method.padics.pAdic_base import pAdicBase
        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
        sage: from modular_method.diophantine_equations.conditions import ConditionalValue
        sage: from modular_method.diophantine_equations.conditions import TreeCondition
        sage: from modular_method.diophantine_equations.conditions import conditional_over_values
        sage: R.<x> = ZZ[]
        sage: C = PolynomialCondition(x^2 + 1)
        sage: pAdics = pAdicBase(QQ, 5)
        sage: def f(val, con):
        ....:     Y, N = C.pAdic_tree(pAdic_tree=con.pAdic_tree(pAdics=pAdics),
        ....:                         complement=True, precision=val)
        ....:     success = "Solution for precision " + str(val)
        ....:     return ConditionalValue([(success, TreeCondition(Y)),
        ....:                              (None, TreeCondition(N))])
        ....: 
        sage: def g(val, con):
        ....:     con = (CoprimeCondition([x]) if con is None else con)
        ....:     return f(val, con)
        ....: 
        sage: conditional_over_values(f, [4, 3, 2, 1])
        Traceback (most recent call last):
        ...
        AttributeError: 'NoneType' object has no attribute 'pAdic_tree'
        sage: conditional_over_values(g, [4, 3, 2, 1])
        Solution for precision 4 if x == 182, 443 mod 625
        Solution for precision 3 if x is 1 of 8 possibilities mod 625
        Solution for precision 2 if x is 1 of 8 possibilities mod 125
        Solution for precision 1 if x == 2, 3, 8, 12, 13, 17, 22, 23 mod 25
        None                     if x == 0, 1, 4 mod 5
    
    The given function can have multiple answers which are combined into one::

        sage: from modular_method.padics.pAdic_base import pAdicBase
        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: from modular_method.diophantine_equations.conditions import PolynomialCondition
        sage: from modular_method.diophantine_equations.conditions import ConditionalValue
        sage: from modular_method.diophantine_equations.conditions import TreeCondition
        sage: from modular_method.diophantine_equations.conditions import conditional_over_values
        sage: R.<x> = ZZ[]
        sage: C1 = PolynomialCondition(x^4 + 4)
        sage: C2 = PolynomialCondition(x^2 + 4)
        sage: pAdics = pAdicBase(QQ, 13)
        sage: def f(val, con):
        ....:     Y1, N1 = C1.pAdic_tree(pAdic_tree=con.pAdic_tree(pAdics=pAdics),
        ....:                            complement=True, precision=val)
        ....:     success1 = "Solution for x^4 + 1 with precision " + str(val)
        ....:     Y2, N2 = C2.pAdic_tree(pAdic_tree=N1, complement=True,
        ....:                            precision=val)
        ....:     success2 = "Solution for x^2 + 1 with precision " + str(val)
        ....:     return ConditionalValue([(success1, TreeCondition(Y1)),
        ....:                              (success2, TreeCondition(Y2)),
        ....:                              (None, TreeCondition(N2))])
        ....: 
        sage: start = CoprimeCondition([x])
        sage: conditional_over_values(f, [2, 1], start_condition=start)
        Solution for x^4 + 1 with precision 2 if x == 69, 71, 98, 100 mod 169
        Solution for x^2 + 1 with precision 2 if x == 29, 140 mod 169
        Solution for x^4 + 1 with precision 1 if x is 1 of 48 possibilities mod 169
        Solution for x^2 + 1 with precision 1 if x is 1 of 24 possibilities mod 169
        None                                  if x == 0, 1, 2, 5, 8, 11, 12 mod 13

    """
    vals = [None]
    cons = [start_condition]
    values = iter(values)
    for value in values:
        if not (None in vals):
            break
        i = vals.index(None)
        con = cons.pop(i)
        vals.pop(i)
        result = function(value, con)
        if not isinstance(result, ConditionalValue):
            result = [(result, con)]
        for fval, fcon in result:
            if not fcon.never():
                try:
                    i = vals.index(fval)
                    cons[i] = cons[i] | fcon
                except ValueError:
                    vals.append(fval)
                    cons.append(fcon)
    if not singleton and len(vals) == 1:
        return vals[0]
    else:
        return ConditionalValue(list(zip(vals, cons)))
        
def conditional_product(*args):
    r"""Create a single ConditionalValue from multiple ConditionalValues.

    INPUT:
    
    Any amount of arguments, all of which can be an instance of
    ConditionalValue or any other value.

    OUTPUT:

    A ConditionalValue of which the values are all possible lists with
    as i-th entry a possible value of the i-th given
    ConditionalValue. For each value the corresponding condition is
    the condition that all the corresponding conditions for each entry
    in the respective ConditionalValue hold.

    Note that if each entry given is not a ConditionalValue, the
    result will also not be a ConditionalValue.

    EXAMPLE::

        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition, PolynomialCondition, ConditionalValue, conditional_product
        sage: R.<x, y> = ZZ[]
        sage: C1 = CoprimeCondition([x, y])
        sage: V1 = ConditionalValue([(2, C1), (7, ~C1)])
        sage: C2 = PolynomialCondition(x^2 + y^2 - 4)
        sage: V2 = ConditionalValue([(27, C2), (0, ~C2)])
        sage: conditional_product(V1, V2)
        (2, 27) if ('x', 'y') are pairwise coprime and x^2 + y^2 - 4 == 0
        (2, 0)  if ('x', 'y') are pairwise coprime and x^2 + y^2 - 4 ~= 0
        (7, 27) if ('x', 'y') are not pairwise coprime and x^2 + y^2 - 4 == 0
        (7, 0)  if ('x', 'y') are not pairwise coprime and x^2 + y^2 - 4 ~= 0

    """
    if len(args) == 0:
        raise ValueError("conditional_product requires at least one argument.")
    args = [(arg if isinstance(arg, ConditionalValue) else [(arg, None)])
            for arg in args]
    # Helper function
    def combine_conditions(C1, C2):
        if C1 is None:
            return C2
        if C2 is None:
            return C1
        return C1 & C2
    # Make a list of lists of value condition pairs
    result = itertools.product(*args)
    # Turn into a list of pairs of a list of values
    # and the corresponding list of conditions.
    result = [zip(*val_con) for val_con in result]
    # Turn into a list of pairs of a list of values
    # and the corresponding condition, combined with and.
    result = [(val, reduce(combine_conditions, con, None))
              for val, con in result]
    if len(result) == 1:
        return result[0][0]
    else:
        return ConditionalValue(result)
