r"""Classes to work with Frey-Hellegouarch curves.

A Frey-Hellegouarch curve, often referred to simply as a Frey curve,
is an elliptic curve of which the coefficients depend on an (unknown)
solution of a diophantine equation. For our purposes a Frey curve will
be an elliptic curve defined over a ring $R$, that is a (multivariate)
polynomial ring over some number field $L$. The variables of $R$ are
assumed to take on undetermined values in the ring of integers of some
subfield $K$ of $L$, which might satisfy some constraints. These
variables are known as the parameters of the Frey curve.

This file provides two variants of Frey curves. The class
:class:`FreyCurve` provides the basic implementation of a Frey
curve. However computing newforms for such a Frey curve only works if
its defined over the rationals. The class :class:`FreyQcurve` extends
the class :class:`FreyCurve` and the class :class:`Qcurve` and
provides a way to work with Frey curves that are also $\Q$-curves.

EXAMPLES:

The classical example, Fermat's Last Theorem that $x^n + y^n = z^n$
has no non-trivial solutions for $n > 3$. In this example we use a
Frey curve based on a solution $(a, b, c)$ with $a b c \ne 0$ for $n$
a prime number $l > 2$. Without loss of generality we assume that $b
\equiv 1$ modulo 4. In this example we use `al`, `bl`, and `cl` for
$a^l$, $b^l$, and $c^l$ respectively::

    sage: from modular_method.diophantine_equations.conditions import PowerCondition
    sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
    sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
    sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
    sage: R.<al, bl> = QQ[]
    sage: cl = al + bl
    sage: C = (PowerCondition(al, 3) & PowerCondition(bl, 3) & PowerCondition(cl, 3) &
    ....:      CongruenceCondition(bl - 1, 4) & CoprimeCondition([al, bl]))
    sage: E = FreyCurve([0, bl - al, 0, -al*bl, 0], condition=C); E
    Frey curve defined by y^2 = x^3 + (-al+bl)*x^2 + (-al*bl)*x over Rational Field with parameters (al, bl)
    sage: E.minimal_model(2)
    Frey curve defined by y^2 = x^3 + (-al+bl)*x^2 + (-al*bl)*x over Rational Field with parameters (al, bl)                                                                                                                                                      if ('al', 'bl') is 1 of 12 possibilities mod 16
    Elliptic Curve defined by y^2 + x*y + (1/4*al-1/4*bl+1/4)*y = x^3 + (1/2*al-1/2*bl+1/2)*x^2 + (1/16*al^2-3/16*al*bl+1/16*bl^2+1/8*al-1/8*bl+1/16)*x + (-1/64*al^2*bl+1/64*al*bl^2-1/64*al*bl) over Multivariate Polynomial Ring in al, bl over Rational Field if ('al', 'bl') is 1 of 4 possibilities mod 16
    sage: E.discriminant().factor()
    (16) * bl^2 * al^2 * (al + bl)^2
    sage: E.conductor()
    Warning: Assuming that al and bl are coprime.
    2^n0*Rad_P( (16) * bl^2 * al^2 * (al + bl)^2 )
     where 
    n0 = 4 if ('al', 'bl') == (3, 5), (7, 1) mod 8
         3 if ('al', 'bl') is 1 of 4 possibilities mod 16
         0 if ('al', 'bl') is 1 of 8 possibilities mod 32
         1 if ('al', 'bl') is 1 of 8 possibilities mod 32
    sage: E.newform_candidates()
    Warning: The bad primes chosen by default only take into account primes of additive reduction.
    []

A more involved example for the equation $x^n + y^n = 2 z^n$ from the
article "Winding quotients and some variants of Fermat's Last Theorem"
by Henri Darmon and Loïc Merel (1997). We let $(a, b, c)$ be a
solution of this equation for $n$ a prime number $p \ge 7$ with
$\gcd(a, b, c) = 1$ and $|a b c| > 1$. We denote by `ap`, `bp`, and
`cp` the values $a^p$, $b^p$ and $c^p$ respectively. As in the article
we assume that $a \equiv -1$ modulo 4::

    sage: from modular_method.diophantine_equations.conditions import PowerCondition
    sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
    sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
    sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
    sage: R.<ap, cp> = QQ[]
    sage: bp = 2*cp - ap
    sage: C = (PowerCondition(ap, 7) & PowerCondition(bp, 7) & PowerCondition(cp, 7) &
    ....:      CongruenceCondition(ap + 1, 4) & CoprimeCondition([ap, cp]))
    sage: E = FreyCurve([0, -ap - 2*cp, 0, 2*ap*cp, 0], condition=C); E
    Frey curve defined by y^2 = x^3 + (-ap-2*cp)*x^2 + 2*ap*cp*x over Rational Field with parameters (ap, cp)
    sage: E.discriminant().factor()
    (64) * cp^2 * ap^2 * (ap - 2*cp)^2
    sage: E.conductor(additive_primes=[2])
    2^n0*Rad_P( (64) * cp^2 * ap^2 * (ap - 2*cp)^2 )
     where 
    n0 = 5 if ('ap', 'cp') == (3, 1), (3, 3) mod 4
         1 if ('ap', 'cp') is 1 of 32 possibilities mod 128
    sage: E.newform_candidates(bad_primes=[2])
    [q - 2*q^5 + O(q^6)] if ('ap', 'cp') == (3, 1), (3, 3) mod 4
    []                   if ('ap', 'cp') is 1 of 32 possibilities mod 128

We use the Frey curve for the equation $A^4 + B^2 = C^p$ from the
article "Galois representations attached to Q-curves and the
generalized Fermat equation $A^4 + B^2 = C^p$" by Jordan S. Ellenberg
(2004) as an example of a Frey Q-curve. We will let $(A, B, C)$ be a
solution to this equation for $p \ge 211$ prime with $ A B \ne 0 $ and
$B \not\equiv 1$ modulo 4. We will use `Cp` to denote $C^p$::

    sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
    sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
    sage: from modular_method.elliptic_curves.frey_curves import FreyQcurve
    sage: R.<A, B> = QQ[]
    sage: Cp = A^4 + B^2
    sage: con = CoprimeCondition([A, B]) & ~CongruenceCondition(B - 1, 4)
    sage: K.<i> = QuadraticField(-1)
    sage: a_invariants = [0, 2*(1+i)*A, 0, B + i*A^2, 0]
    sage: E = FreyQcurve(a_invariants, condition=con, guessed_degrees=[2])
    sage: E.discriminant().factor()
    ((-64*i)) * (A^2 + (i)*B) * (A^2 + (-i)*B)^2
    sage: E.decomposition_field()
    Number Field in i with defining polynomial x^2 + 1 with i = 1*I
    sage: E.does_decompose()
    True
    sage: E.newform_candidates(bad_primes=K.primes_above(2))
    [q - 2*q^3 + O(q^6), q - 4*q^5 + O(q^6), q + 4*q^5 + O(q^6), q + 2*q^3 + O(q^6), q + 1/2*a4*q^3 + O(q^6), q + 1/2*a4*q^3 + O(q^6)] if ('A', 'B') is 1 of 6 possibilities mod 4
    [q - 2*q^5 + O(q^6)]                                                                                                               if ('A', 'B') == (0, 3), (2, 3) mod 4

AUTHORS:

- Joey van Langen (2019-03-06): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy

from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing
from sage.rings.morphism import RingHomomorphism_from_base

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
from sage.schemes.elliptic_curves.constructor import EllipticCurve

from sage.rings.number_field.number_field import is_NumberField

from sage.all import ZZ, QQ, Integer
from sage.misc.cachefunc import cached_method
from sage.arith.functions import lcm
from sage.arith.misc import gcd

from sage.misc.misc_c import prod as product

from modular_method.elliptic_curves.Qcurves import Qcurve_base as Qcurve
from modular_method.elliptic_curves.Qcurves import _rational_maps_of_urst
from modular_method.elliptic_curves.Qcurves import _scalar_of_rational_maps
from modular_method.elliptic_curves.tates_algorithm import tates_algorithm
from modular_method.elliptic_curves.twist import is_twist

from modular_method.padics.pAdic_base import pAdicBase
from modular_method.padics.pAdic_tree import pAdicTree, pAdicNode

from modular_method.number_fields.field_constructors import _write_as_im_gen_map
from modular_method.number_fields.field_constructors import _concat_maps
from modular_method.number_fields.field_constructors import write_as_extension
from modular_method.number_fields.field_constructors import composite_field
from modular_method.number_fields.galois_group import galois_field_change

from modular_method.diophantine_equations.conditions import TextCondition
from modular_method.diophantine_equations.conditions import TreeCondition
from modular_method.diophantine_equations.conditions import ConditionalValue
from modular_method.diophantine_equations.conditions import ConditionalExpression
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_over_values

from modular_method.modular_forms.newform_wrapper import get_newforms

class FreyCurve(EllipticCurve_generic):
    r"""A Frey-Hellegouarch curve.

    A Frey-Hellegouarch curve, or simply a Frey curve, is an elliptic
    curve defined over a (multivariate) polynomial ring $R$ defined
    over some number field $L$. The variables of $R$ are assumed to
    have undetermined values in the ring of integers of some subfield
    $K$ of $L$ and are considered to be parameters of $E$. Therefore
    this elliptic curve is considered as an elliptic curve over $L$.

    This class provides some functionality to compute with this
    elliptic curve as if it is defined over the number field
    $L$. Since the answer of some of these computations might depend
    on the values of the parameters, the result is often given as
    :class:`ConditionalValue`. Most methods also given te option of
    providing additional constraints on the parameters in the form of
    instances of :class:`Condition` and this curve can also internally
    stores a condition on the parameters which is used as a default
    restraint.

    EXAMPLES:

    The classical example, Fermat's Last Theorem that $x^n + y^n = z^n$
    has no non-trivial solutions for $n > 3$. In this example we use a
    Frey curve based on a solution $(a, b, c)$ with $a b c \ne 0$ for $n$
    a prime number $l > 2$. Without loss of generality we assume that $b
    \equiv 1$ modulo 4. In this example we use `al`, `bl`, and `cl` for
    $a^l$, $b^l$, and $c^l$ respectively::

        sage: from modular_method.diophantine_equations.conditions import PowerCondition
        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
        sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
        sage: R.<al, bl> = QQ[]
        sage: cl = al + bl
        sage: C = (PowerCondition(al, 3) & PowerCondition(bl, 3) & PowerCondition(cl, 3) &
        ....:      CongruenceCondition(bl - 1, 4) & CoprimeCondition([al, bl]))
        sage: E = FreyCurve([0, bl - al, 0, -al*bl, 0], condition=C); E
        Frey curve defined by y^2 = x^3 + (-al+bl)*x^2 + (-al*bl)*x over Rational Field with parameters (al, bl)
        sage: E.minimal_model(2)
        Frey curve defined by y^2 = x^3 + (-al+bl)*x^2 + (-al*bl)*x over Rational Field with parameters (al, bl)                                                                                                                                                      if ('al', 'bl') is 1 of 12 possibilities mod 16
        Elliptic Curve defined by y^2 + x*y + (1/4*al-1/4*bl+1/4)*y = x^3 + (1/2*al-1/2*bl+1/2)*x^2 + (1/16*al^2-3/16*al*bl+1/16*bl^2+1/8*al-1/8*bl+1/16)*x + (-1/64*al^2*bl+1/64*al*bl^2-1/64*al*bl) over Multivariate Polynomial Ring in al, bl over Rational Field if ('al', 'bl') is 1 of 4 possibilities mod 16
        sage: E.discriminant().factor()
        (16) * bl^2 * al^2 * (al + bl)^2
        sage: E.conductor()
        Warning: Assuming that al and bl are coprime.
        2^n0*Rad_P( (16) * bl^2 * al^2 * (al + bl)^2 )
         where 
        n0 = 4 if ('al', 'bl') == (3, 5), (7, 1) mod 8
             3 if ('al', 'bl') is 1 of 4 possibilities mod 16
             0 if ('al', 'bl') is 1 of 8 possibilities mod 32
             1 if ('al', 'bl') is 1 of 8 possibilities mod 32
        sage: E.newform_candidates()
        Warning: The bad primes chosen by default only take into account primes of additive reduction.
        []

    A more involved example for the equation $x^n + y^n = 2 z^n$ from the
    article "Winding quotients and some variants of Fermat's Last Theorem"
    by Henri Darmon and Loïc Merel (1997). We let $(a, b, c)$ be a
    solution of this equation for $n$ a prime number $p \ge 7$ with
    $\gcd(a, b, c) = 1$ and $|a b c| > 1$. We denote by `ap`, `bp`, and
    `cp` the values $a^p$, $b^p$ and $c^p$ respectively. As in the article
    we assume that $a \equiv -1$ modulo 4::

        sage: from modular_method.diophantine_equations.conditions import PowerCondition
        sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
        sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
        sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
        sage: R.<ap, cp> = QQ[]
        sage: bp = 2*cp - ap
        sage: C = (PowerCondition(ap, 7) & PowerCondition(bp, 7) & PowerCondition(cp, 7) &
        ....:      CongruenceCondition(ap + 1, 4) & CoprimeCondition([ap, cp]))
        sage: E = FreyCurve([0, -ap - 2*cp, 0, 2*ap*cp, 0], condition=C); E
        Frey curve defined by y^2 = x^3 + (-ap-2*cp)*x^2 + 2*ap*cp*x over Rational Field with parameters (ap, cp)
        sage: E.discriminant().factor()
        (64) * cp^2 * ap^2 * (ap - 2*cp)^2
        sage: E.conductor(additive_primes=[2])
        2^n0*Rad_P( (64) * cp^2 * ap^2 * (ap - 2*cp)^2 )
         where 
        n0 = 5 if ('ap', 'cp') == (3, 1), (3, 3) mod 4
             1 if ('ap', 'cp') is 1 of 32 possibilities mod 128
        sage: E.newform_candidates(bad_primes=[2])
        [q - 2*q^5 + O(q^6)] if ('ap', 'cp') == (3, 1), (3, 3) mod 4
        []                   if ('ap', 'cp') is 1 of 32 possibilities mod 128

    """
    def __init__(self, curve, parameter_ring=ZZ, condition=None):
        r"""Initialize a Frey curve

        INPUT:

        - ``curve`` -- An elliptic curve or any argument that would
          produce such a curve when passed to the constructor
          :func:`EllipticCurve`. The elliptic curve should be defined
          over some (multivariate) polynomial ring $R$ which in turn
          is defined over some number field $L$. This will be the Frey
          curve.

        - ``parameter_ring`` -- The ring of integers of a subfield $K$
          of $L$ (default: ZZ). This is the ring in which the
          variables of the polynomial ring over which this curve is
          defined can take values.

        - ``condition`` -- An instance of :class:`Condition` or None
          (default: None) giving a condition which must hold for the
          values of the variables of $R$. If set to None will assume
          that all values for these variables are allowed.

        EXAMPLES:

        The classical example, Fermat's Last Theorem that $x^n + y^n = z^n$
        has no non-trivial solutions for $n > 3$. In this example we use a
        Frey curve based on a solution $(a, b, c)$ with $a b c \ne 0$ for $n$
        a prime number $l > 2$. Without loss of generality we assume that $b
        \equiv 1$ modulo 4. In this example we use `al`, `bl`, and `cl` for
        $a^l$, $b^l$, and $c^l$ respectively::

            sage: from modular_method.diophantine_equations.conditions import PowerCondition
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
            sage: R.<al, bl> = QQ[]
            sage: cl = al + bl
            sage: C = (PowerCondition(al, 3) & PowerCondition(bl, 3) & PowerCondition(cl, 3) &
            ....:      CongruenceCondition(bl - 1, 4) & CoprimeCondition([al, bl]))
            sage: FreyCurve([0, bl - al, 0, -al*bl, 0], condition=C)
            Frey curve defined by y^2 = x^3 + (-al+bl)*x^2 + (-al*bl)*x over Rational Field with parameters (al, bl)

        A more involved example for the equation $x^n + y^n = 2 z^n$ from the
        article "Winding quotients and some variants of Fermat's Last Theorem"
        by Henri Darmon and Loïc Merel (1997). We let $(a, b, c)$ be a
        solution of this equation for $n$ a prime number $p \ge 7$ with
        $\gcd(a, b, c) = 1$ and $|a b c| > 1$. We denote by `ap`, `bp`, and
        `cp` the values $a^p$, $b^p$ and $c^p$ respectively. As in the article
        we assume that $a \equiv -1$ modulo 4::

            sage: from modular_method.diophantine_equations.conditions import PowerCondition
            sage: from modular_method.diophantine_equations.conditions import CoprimeCondition
            sage: from modular_method.diophantine_equations.conditions import CongruenceCondition
            sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
            sage: R.<ap, cp> = QQ[]
            sage: bp = 2*cp - ap
            sage: C = (PowerCondition(ap, 7) & PowerCondition(bp, 7) & PowerCondition(cp, 7) &
            ....:      CongruenceCondition(ap + 1, 4) & CoprimeCondition([ap, cp]))
            sage: FreyCurve([0, -ap - 2*cp, 0, 2*ap*cp, 0], condition=C)
            Frey curve defined by y^2 = x^3 + (-ap-2*cp)*x^2 + 2*ap*cp*x over Rational Field with parameters (ap, cp)

        It is not required to give a condition::

            sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
            sage: R.<al, bl> = QQ[]
            sage: cl = al + bl
            sage: FreyCurve([0, bl - al, 0, -al*bl, 0])
            Frey curve defined by y^2 = x^3 + (-al+bl)*x^2 + (-al*bl)*x over Rational Field with parameters (al, bl)

        Note that the elliptic curve does not necessarily have to be
        associated to a diophantine equation::

            sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
            sage: R.<a, b> = QQ[]
            sage: FreyCurve([a-1, b + 2, a^3 + a*b, b^2 + 4*b, b^3])
            Frey curve defined by y^2 + (a-1)*x*y + (a^3+a*b)*y = x^3 + (b+2)*x^2 + (b^2+4*b)*x + b^3 over Rational Field with parameters (a, b)

        You can construct Frey curves over any number field as long as
        the parameter ring has a natural inclusion into that field::

            sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
            sage: R.<A, B> = QQ[]
            sage: Cp = A^4 + B^2
            sage: con = CoprimeCondition([A, B]) & ~CongruenceCondition(B - 1, 4)
            sage: K.<i> = QuadraticField(-1)
            sage: a_invariants = [0, 2*(1+i)*A, 0, B + i*A^2, 0]
            sage: FreyCurve(a_invariants, condition=con)
            Frey curve defined by y^2 = x^3 + ((2*i+2)*A)*x^2 + ((i)*A^2+B)*x over Number Field in i with defining polynomial x^2 + 1 with i = 1*I with parameters (A, B)
            sage: FreyCurve(a_invariants, parameter_ring=K.ring_of_integers(), condition=con)
            Frey curve defined by y^2 = x^3 + ((2*i+2)*A)*x^2 + ((i)*A^2+B)*x over Number Field in i with defining polynomial x^2 + 1 with i = 1*I with parameters (A, B)

        """
        if not isinstance(curve, EllipticCurve_generic):
            curve = EllipticCurve(curve)
        S = curve.base_ring()
        if not (is_PolynomialRing(S) or is_MPolynomialRing(S)):
            raise ValueError("The coefficient ring " + str(S) +
                             " is not a polynomial ring.")
        base = S.base_ring()
        EllipticCurve_generic.__init__(self, S, curve.a_invariants())
        self._R = parameter_ring
        if self._R == ZZ:
            self._R_to_base = QQ.embeddings(base)[0]
        else:
            self._R_to_base = self._R.number_field().embeddings(base)[0]
        self._condition = condition

    def definition_field(self):
        r"""Give the field over which this Frey curve is defined.

        Even though the Frey curve is defined over some (multivariate)
        polynomial ring $R$ over some number field $L$, since the
        variables of $R$ are assumed to have values in some subfield
        $K$ of $L$ the curve can be assumed to be defined over $L$.

        OUTPUT:
        
        The base ring of the polynomial ring over which this Frey
        curve is defined.

        .. SEE_ALSO::

            :meth:`base_ring`,
            :meth:`parameters`

        """
        return self.base_ring().base()
    
    def parameters(self):
        r"""Give the parameters on which this Frey curve depends.

        OUTPUT:

        The variables of the polynomial ring over which this Frey
        curve is defined.

        .. SEE_ALSO::

            :meth:`base_ring`,
            :meth:`definition_field`

        """
        return self.base().gens()

    def parameter_ring(self):
        r"""Give the ring in which the parameters can take values.

        OUTPUT:

        A ring in which the parameters of this Frey curve take
        values.

        """
        return self._R

    @cached_method
    def primes_of_possible_additive_reduction(self):
        r"""Compute the primes at which this curve could have additive
        reduction.

        Tries to find all the primes that divide both the discriminant
        $D$ and the invariant $c_4$ of this Frey curve. Since both of
        these invariants are elements of a (multivariate) polynomial
        ring $R$ these primes are not necessarily uniquely determined.
        In case an additional condition is needed on the values of the
        variables of the polynomials, a warning will be printed with
        the additional condition required.

        OUTPUT:

        A list of primes of the definition field of this Frey
        curve. In case the definition field is $\QQ$ these will be
        prime numbers. In case it is a number field they will be
        maximal ideals of its ring of integers. In any case they
        satisfy one of the following conditions.

        - If the curve has no parameters, the list will contain all the primes
          dividing the greatest common divisor of $D$ and $c_4$

        - If the curve has one parameter, the list will contain all the
          primes dividing the resultant of $D$ and $c_4$

        - If the curve has two parameters and $c4$ and $D$ are
          homogeneous in the parameters, the parameters will be
          assumed to be coprime and the list will contain all the
          primes dividing the Macaulay resultant of $D$ and $c_4$.
        
        - If the curve has more than two parameters, the list will
          contain all the primes dividing all the coefficients of $D$
          and $c_4$. Furthermore $D$ and $c_4$ will be assumed to be
          coprime.

        """
        R = self.base()
        n = len(R.gens())
        K = R.base()
        c4 = self.c4()
        D = self.discriminant()
        if K.is_subring(QQ):
            if n == 1:
                return QQ(c4.resultant(D)).numerator().prime_factors()
            elif n == 2 and c4.is_homogeneous() and D.is_homogeneous():
                print("Warning: Assuming that " + str(self.parameters()[0]) +
                       " and " + str(self.parameters()[1]) + " are coprime.")
                return QQ(c4.macaulay_resultant(D)).numerator().prime_factors()
            else:
                N = lcm(lcm(QQ(c).denominator() for c in c4.coefficients()),
                        lcm(QQ(c).denominator() for c in D.coefficients()))
                M = gcd(gcd([ZZ(c) for c in (N * c4).coefficients()]),
                        gcd([ZZ(c) for c in (N * D).coefficients()]))
                result = QQ(M/N).numerator().prime_factors()
                print("Warning: Assuming that %s and %s "%(c4,D) +
                       " are coprime outside %s."%(tuple(result),))
                return result
        if n == 1:
            return K.ideal(c4.resultant(D)).prime_factors()
        elif n == 2 and c4.is_homogeneous() and D.is_homogeneous():
            print("Warning: Assuming that " + str(self.parameters()[0]) +
                   " and " + str(self.parameters()[1]) + " are coprime.")
            return K.ideal(c4.macaulay_resultant(D)).prime_factors()
        else:
            I = sum(K.ideal(c) for c in c4.coefficients())
            J = sum(K.ideal(c) for c in D.coefficients())
            result = (I + J).prime_factors()
            print("Warning: Assuming that %s and %s"%(c4,D) +
                   " are coprime outside %s."%(tuple(P._repr_short()
                                                    for P in result),))
            return result

    @cached_method(key=(lambda self, p, c, v, pc:
                        (p, (self._condition if c is None else c), pc)))
    def _initial_tree(self, prime, condition=None, verbose=False, precision_cap=20):
        r"""Give a p-adic tree of possible values for the parameters.

        INPUT:

        - ``prime`` -- A (maximal) prime ideal of the ring in which
          the parameters take values or any generator thereof if it is
          principal.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.
        
        OUTPUT:

        A p-adic tree containing all the possible p-adic values of the
        the parameters satisfying the given condition.

        """
        Tfull = pAdicTree(variables=self.parameters(),
                          pAdics=pAdicBase(self._R, prime))
        if condition is None:
            condition = self._condition
        if condition is None:
            return Tfull
        else:
            return condition.pAdic_tree(pAdic_tree=Tfull, verbose=verbose,
                                        precision_cap=precision_cap)

    @cached_method(key=(lambda self, p, c, v, pc:
                        (p, (self._condition if c is None else c), pc)))
    def local_data(self, prime, condition=None, verbose=False,
                   precision_cap=20):
        r"""Give the local data of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:
        
        The local data of this Frey curve at the given prime under the
        given condition on the parameters. This could be a conditional
        value as the local data might depend on the value of the
        parameters in this Frey curve.

        .. SEEALSO::

            :meth:`local_data`,
            :meth:`kodaira_symbol`,
            :meth:`minimal_model`,
            :meth:`conductor_exponent`,
            :meth:`reduction_type`,

        """
        pAdics = pAdicBase(self.definition_field(), prime)
        Tp = self._initial_tree(pAdics.prime_below(self._R),
                                condition=condition,
                                verbose=(verbose-1 if verbose>0 else verbose))
        result = tates_algorithm(self, initial_values=Tp,
                                 coefficient_ring=self.base(), pAdics=pAdics,
                                 verbose=verbose, precision_cap=precision_cap)
        if len(result) == 1:
            return result[0][0]
        else:
            return result

    @cached_method(key=(lambda self, p, iso, c, v, pc:
                        (p, iso, (self._condition if c is None else c), pc)))
    def minimal_model(self, prime, isomorphism=False, condition=None,
                      verbose=False, precision_cap=20):
        r"""Give a minimal model of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``isomorphism`` -- A boolean value (default: False). If set
          to True this method will also return the change of
          weierstrass model necessary to change this curve into its
          minimal model and the change needed to change the minimal
          model into this curve.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:
        
        An elliptic curve that is a minimal model of this curve at the
        given prime under the given condition on the parameters.

        If `isomorphism` is set to True then the return value will be
        a tuple containing such a minimal model as a first value. The
        second value will be a tuple consisting of two changes of
        weierstrass models, the first changing this curve into the
        minimal model and the second changing the minimal model back
        into this curve. Each of these change of variables is given as
        a tuple (u, r, s, t), where u, r, s, and t are the values that
        need to be passed to :meth:`change_weierstrass_model` to enact
        the correct change.

        If the return value might differ depending on several
        sub conditions of the given condition, then this method will
        return a ConditionalValue consisting of the corresponding
        output for each of these conditions.

        .. SEEALSO::

            :meth:`local_data`,

        """
        if self.local_data.is_in_cache(prime, condition, verbose,
                                       precision_cap):
            local_data = self.local_data(prime, condition, verbose,
                                         precision_cap)
            if isomorphism:
                return apply_to_conditional_value(lambda x:
                                                  (x.minimal_model(),
                                                   x._urst,
                                                   x._urst_inv),
                                                  local_data)
            else:
                return apply_to_conditional_value(lambda x: x.minimal_model(),
                                                  local_data)
        pAdics = pAdicBase(self.definition_field(), prime)
        Tp = self._initial_tree(pAdics.prime_below(self._R), condition,
                                verbose=(verbose-1 if verbose>0 else verbose))
        calc = ['minimal_model']
        if isomorphism:
            calc.append('isomorphism')
        result = tates_algorithm(self, initial_values=Tp,
                                 coefficient_ring=self.base(), pAdics=pAdics,
                                 verbose=verbose, precision_cap=precision_cap,
                                 only_calculate=calc)
        if isomorphism:
            return apply_to_conditional_value(lambda x: tuple(x),
                                              result)
        else:
            return apply_to_conditional_value(lambda x: x[0],
                                              result)

    def reduction(self, prime, model=None, condition=None,
                  verbose=False, precision_cap=20):
        r"""Give the reduction of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``model`` -- A Weierstrass model of this elliptic curve or
          None (default) giving the model of which the reduction
          should be computed. If set to None will use a minimal model
          of this curve at the given prime instead.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:
        
        An elliptic curve that is the reduction of the given model of
        this curve at a given prime.

        If the return value depends on the specific value of the
        parameters, will instead return a ConditionalValue containing
        the different reductions with the corresponding conditions on
        the parameters for which they are the right reduction.

        If any value of the parameters satisfying the given condition
        would not give an elliptic curve as a reduction, this function
        will produce an error.

        .. SEEALSO::

            :meth:`minimal_model`,

        """
        if condition == None:
            condititon = self._condition
        if model == None:
            model = self.minimal_model(prime, condition=condition,
                                       verbose=(verbose
                                                if verbose <= 0
                                                else verbose - 1),
                                       precision_cap=precision_cap)
        if isinstance(model, ConditionalValue):
            return apply_to_conditional_value(lambda E, C:
                                              self.reduction(prime,
                                                             model=E, condition=C,
                                                             verbose=verbose,
                                                             precision_cap=precision_cap),
                                              model,
                                              use_condition=True,
                                              default_condition=condition)
        pAdics = pAdicBase(self.definition_field(), prime)
        pAdicsR = pAdics.pAdics_below(self._R)
        T = C.pAdic_tree(pAdics=pAdicsR,
                         verbose=(verbose if verbose <= 0 else verbose - 1),
                         precision_cap=precision_cap).root()
        ext = pAdics.extension_multiplicity(pAdicsR)
        E_val = ZZ(min(pAdics.valuation(cf) for a in model.a_invariants()
                       for cf in a.coefficients())) // ext
        T_val = T.minimum_full_level()
        val = max(1, min(1 - E_val, T_val, precision_cap))
        F = pAdics.residue_field()
        result = {}
        for N in T.children_at_level(val):
            vals = N.representative()
            E = EllipticCurve([F(a(vals))
                               for a in model.a_invariants()])
            if not (E in result):
                result[E] = pAdicNode(pAdics=T.pAdics(), width=T.width)
            result[E].merge(N)
        result = [(val, pAdicTree(self.parameters(), root=T))
                  for T, val in result.items()]
        return ConditionalValue([(val, TreeCondition(T))
                                 for val, T in result])

    @cached_method(key=(lambda self, p, c, v, pc:
                        (p, (self._condition if c is None else c), pc)))
    def kodaira_symbol(self, prime, condition=None, verbose=False,
                       precision_cap=20):
        r"""Give the kodaira symbol of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:
        
        The KodairaSymbol representing the type of reduction of this
        curve at the given prime for the parameters satisfying the
        given condition. This could be a conditional value as it might
        depend on the value of the parameters in this curve.

        .. SEEALSO::

            :meth:`local_data`,

        """
        if self.local_data.is_in_cache(prime, condition, verbose,
                                       precision_cap):
            local_data = self.local_data(prime, condition, verbose,
                                         precision_cap)
            return apply_to_conditional_value(lambda x: x.kodaira_symbol(),
                                              local_data)
        
        pAdics = pAdicBase(self.definition_field(), prime)
        Tp = self._initial_tree(pAdics.prime_below(self._R),
                                condition=condition,
                                verbose=(verbose-1 if verbose>0 else verbose))
        result = tates_algorithm(self, initial_values=Tp,
                                 coefficient_ring=self.base(), pAdics=pAdics,
                                 verbose=verbose, precision_cap=precision_cap,
                                 only_calculate=['minimal_model'])
        if len(result) == 1:
            return result[0][0][0]
        else:
            return ConditionalValue([(val[0], con) for val, con in result])
    
    @cached_method(key=(lambda self, p, c, v, pc:
                        (p, (self._condition if c is None else c), pc)))
    def conductor_exponent(self, prime, condition=None, verbose=False,
                           precision_cap=20):
        r"""Give the conductor exponent of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        The exponent of the conductor of this Frey curve at the given
        prime for parameters satisfying the given condition. This
        could be a conditional value as the conductor exponent might
        depend on the value of the parameters in this Frey curve.

        .. SEEALSO::

            :meth:`local_data`,
            :meth:`conductor`

        """
        if self.local_data.is_in_cache(prime, condition, verbose,
                                       precision_cap):
            local_data = self.local_data(prime, condition, verbose,
                                         precision_cap)
            return apply_to_conditional_value(lambda x:
                                              x.conductor_valuation(),
                                              local_data)
        
        pAdics = pAdicBase(self.definition_field(), prime)
        Tp = self._initial_tree(pAdics.prime_below(self._R),
                                condition=condition,
                                verbose=(verbose-1 if verbose>0 else verbose))
        result = tates_algorithm(self, initial_values=Tp,
                                 coefficient_ring=self.base(), pAdics=pAdics,
                                 verbose=verbose, precision_cap=precision_cap,
                                 only_calculate=['conductor'])
        if(len(result) == 1):
            return result[0][0][0]
        else:
            return ConditionalValue([(val[0], con) for val,con in result])

    @cached_method(key=(lambda self, p, c, v, pc:
                        (p, (self._condition if c is None else c), pc)))
    def reduction_type(self, prime, condition=None, verbose=False,
                       precision_cap=20):
        r"""Give the reduction type of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        One of the following values

        - None, if the elliptic curve has good reduction at the given
          prime.

        - 1, if the elliptic curve has split multiplicative reduction
          at the given prime.

        - -1, if the elliptic curve has non-split multiplicative
          reduction at the given prime.

        - 0, if the elliptic curve has additive reduction at the given
          prime.

        If there is multiple options for the possible reductions, will
        return a conditional value instead containing the above values
        and the cases in which these occur.

        .. SEEALSO::

            :meth:`local_data`,
            :meth:`has_good_reduction`,
            :meth:`has_bad_reduction`,
            :meth:`has_additive_reduction`,
            :meth:`has_multiplicative_reduction`,
            :meth:`has_split_multiplicative_reduction`,
            :meth:`has_non_split_multiplicative_reduction`

        """
        if self.local_data.is_in_cache(prime, condition, verbose,
                                       precision_cap):
            local_data = self.local_data(prime, condition, verbose,
                                         precision_cap)
            return apply_to_conditional_value(lambda x: x._reduction_type,
                                              local_data)
        pAdics = pAdicBase(self.definition_field(), prime)
        Tp = self._initial_tree(pAdics.prime_below(self._R),
                                condition=condition,
                                verbose=(verbose-1 if verbose>0 else verbose),
                                precision_cap=precision_cap)
        result = tates_algorithm(self, initial_values=Tp,
                                 coefficient_ring=self.base(), pAdics=pAdics,
                                 verbose=verbose, precision_cap=precision_cap,
                                 only_calculate=['reduction_type'])
        if(len(result) == 1):
            return result[0][0][0]
        else:
            return ConditionalValue([(val[0], con) for val,con in result])

    def has_good_reduction(self, prime, condition=None, verbose=False,
                           precision_cap=20):
        r"""Tell whether this curve has good reduction at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        True, if this Frey curve has good reduction at the given
        prime.

        False, if this Frey curve has bad reduction at the given
        prime.

        If the answer depends on the chosen parameters for this Frey
        curve, will return a conditional value containing the above
        values and the conditions for which they occur.

        .. SEEALSO::

            :meth:`reduction_type`

        """
        red_type = self.reduction_type(prime, verbose=verbose,
                                       condition=condition,
                                       precision_cap=precision_cap)
        return apply_to_conditional_value(lambda x: (x == None), red_type)

    def has_bad_reduction(self, prime, condition=None, verbose=False,
                          precision_cap=20):
        r"""Tell whether this curve has bad reduction at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        True, if this Frey curve has bad reduction at the given prime.

        False, if this Frey curve has good reduction at the given
        prime.

        If the answer depends on the chosen parameters for this Frey
        curve will return a conditional value containing the above
        values and the conditions for which they occur.

        .. SEEALSO::

            :meth:`reduction_type`

        """
        red_type = self.reduction_type(prime, verbose=verbose,
                                       condition=condition,
                                       precision_cap=precision_cap)
        return apply_to_conditional_value(lambda x: (x != None), red_type)

    def has_additive_reduction(self, prime, condition=None, verbose=False,
                               precision_cap=20):
        r"""Tell whether this curve has additive reduction at a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        True, if this Frey curve has additive reduction at the given
        prime.

        False, if this Frey curve does not have additive reduction at
        the given prime.

        If the answer depends on the chosen parameters for this Frey
        curve will return a conditional value containing the above
        values and the conditions for which they occur.

        .. SEEALSO::

            :meth:`reduction_type`

        """
        red_type = self.reduction_type(prime, verbose=verbose,
                                       condition=condition,
                                       precision_cap=precision_cap)
        return apply_to_conditional_value(lambda x: (x == 0), red_type)

    def has_split_multiplicative_reduction(self, prime, condition=None,
                                           verbose=False, precision_cap=20):
        r"""Tell whether this curve has split multiplicative reduction at a
        given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        True, if this Frey curve has split multiplicative reduction at
        the given prime.

        False, if this Frey curve does not have split multiplicative
        reduction at the given prime.

        If the answer depends on the chosen parameters for this Frey
        curve will return a conditional value containing the above
        values and the conditions for which they occur.

        .. SEEALSO::

            :meth:`reduction_type`

        """
        red_type = self.reduction_type(prime, verbose=verbose,
                                       condition=condition,
                                       precision_cap=precision_cap)
        return apply_to_conditional_value(lambda x: (x == 1), red_type)

    def has_non_split_multiplicative_reduction(self, prime, condition=None,
                                               verbose=False,
                                               precision_cap=20):
        r"""Tell whether this curve has non-split multiplicative reduction at
        a given prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        True, if this Frey curve has non-split multiplicative
        reduction at the given prime.

        False, if this Frey curve does not have non-split
        multiplicative reduction at the given prime.

        If the answer depends on the chosen parameters for this Frey
        curve will return a conditional value containing the above
        values and the conditions for which they occur.

        .. SEEALSO::

            :meth:`reduction_type`

        """
        red_type = self.reduction_type(prime, verbose=verbose,
                                       condition=condition,
                                       precision_cap=precision_cap)
        return apply_to_conditional_value(lambda x: (x == -1), red_type)

    def has_multiplicative_reduction(self, prime, condition=None,
                                     verbose=False, precision_cap=20):
        r"""Tell whether this curve has multiplicative reduction at a given
        prime.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This should be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the ring of integers
          otherwise.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        True, if this Frey curve has multiplicative reduction at the
        given prime.

        False, if this Frey curve does not have multiplicative
        reduction at the given prime.

        If the answer depends on the chosen parameters for this Frey
        curve will return a ConditionalValue containing the above
        values and the conditions for which they occur.

        .. SEEALSO::

            :meth:`reduction_type`

        """
        red_type = self.reduction_type(prime, verbose=verbose,
                                       condition=condition,
                                       precision_cap=precision_cap)
        return apply_to_conditional_value(lambda x: (x == 1 or x == -1),
                                          red_type)
        
    def base_extend(self, R):
        if (hasattr(R, 'domain') and R.domain() == self.definition_field()):
            dom = self.base_ring()
            codom = dom.change_ring(R.codomain())
            F = RingHomomorphism_from_base(dom.Hom(codom), R)
            result = EllipticCurve_generic.base_extend(self, F)
        else:
            result = EllipticCurve_generic.base_extend(self, R)
        if ((is_PolynomialRing(result.base_ring()) or
             is_MPolynomialRing(result.base_ring())) and
            (result.base_ring().variable_names() ==
             tuple(str(v) for v in self.parameters()))):
            return FreyCurve(result, parameter_ring=self._R,
                             condition=self._condition)
        return result

    def specialize(self, values):
        r"""Give this curve for a specific value of the parameters.

        INPUT:

        - ``values`` -- A list or tuple containing as the i-th entry
          the value that the i-th parameter in should have. The order
          of the parameters is the same as that returned by
          :meth:`parameters`. All values should be elements of the
          ring in which the parameters take value.

        OUTPUT:

        The elliptic curve wherein each parameter is replaced with the
        corresponding value of the given list values.

        .. SEEALSO::

            :meth:`parameter_ring`

        """
        a_invs = [a(tuple(self._R_to_base(val) for val in values))
                  for a in self.a_invariants()]
        return EllipticCurve(a_invs)

    def conductor(self, additive_primes=None, condition=None, verbose=False,
                  precision_cap=20):
        r"""Compute the conductor of this Frey curve.

        INPUT:

        - ``additive_primes`` -- A list containing primes of the
          definition field of this curve or None (default: None). The
          primes in this list should be given as prime number if the
          definition field is $\QQ$ or as maximal ideals
          otherwise. This list should include all the primes at which
          this curve could have additive reduction. If set to None
          will compute this by using the method
          :meth:`primes_of_possible_additive_reduction`.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings.  A higher value will cause
          more messages to be printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        A conditional expression that gives the conductor of this Frey
        curve for each possible value of the parameters on which it
        depends. The left side of this expression is some expression
        that gives the part of the conductor at all the primes given
        in `additive_primes`, whilst the right side is a string
        describing how to compute the part of the conductor that is
        coprime to those primes. The latter contains the operator
        Rad_P which refers to taking the radical of an expression
        ignoring those primes in `additive_primes`.

        """
        if additive_primes is None:
            additive_primes = self.primes_of_possible_additive_reduction()
        result = 1
        for P in additive_primes:
            factor = P^self.conductor_exponent(P, condition=condition,
                                               verbose=verbose,
                                               precision_cap=precision_cap)
            if result == 1:
                result = factor
            elif isinstance(factor, ConditionalExpression):
                operator = ConditionalExpression.PRODUCT_OPERATOR
                result = ConditionalExpression(operator, result, factor)
            else:
                result = result * factor
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR,
                                     result,
                                     "Rad_P( " +
                                     str(self.discriminant().factor()) + " )")
    
    def _trace_of_frobenius(self, pAdics, red_type, condition, verbose,
                            precision_cap):
        r"""Implementation of :meth:`trace_of_frobenius`

        INPUT:

        - `pAdics` -- The p-adics to be used for the Frobenius
          element

        - `red_type` -- The reduction type of this curve

        - `condition` -- A condition on the parameters of this curve

        - `verbose` -- Verbosity argument

        - `precision_cap` -- Bound on the precision of the parameters

        """
        T = condition.pAdic_tree(pAdics=pAdics.pAdics_below(self._R),
                                 verbose=(verbose-1 if verbose>0 else verbose),
                                 precision_cap=precision_cap).root()
        Fp = len(pAdics.residue_field())
        if red_type is None:
            result = {}
            computed = {}
            for N in T.children_at_level(1):
                E = self.specialize(N.representative())
                Ered = E.reduction(pAdics.prime())
                for E2 in computed:
                    a = is_twist(Ered, E2)
                    if a != 0:
                        ap = a * computed[E2]
                        break
                else:
                    Ep = Ered.count_points()
                    ap = 1 + Fp - Ep
                    computed[Ered] = ap
                if ap not in result:
                    result[ap] = pAdicNode(pAdics=N.pAdics(),
                                           width=N.width)
                result[ap].merge(N)
            return result
        elif red_type == 1 or red_type == -1:
            return {red_type*(len(pAdics.residue_field()) + 1): T}
        else:
            raise ValueError("Can not compute trace of frobenius " +
                             "if the curve has additive reduction.")

    def trace_of_frobenius(self, prime, power=1, condition=None,
                           verbose=False, precision_cap=20):
        r"""Compute the trace of a Frobenius element acting on this curve.

        If the elliptic curve has good reduction at the given prime,
        for every prime number $l$ not divisible by that prime the
        $l$-adic galois representation of this curve is unramified at
        that prime and the trace of the Frobenius element at that
        prime is given by this function.

        If the elliptic curve has multiplicative reduction at the
        given prime, for every prime number $l$ not divisible by the
        prime, the mod $l$ galois representation of this curve is
        unramified at that prime if and only if the valuation of the
        discriminant is divisible by $l$. This function returns the
        trace of the Frobenius element at that prime in such a case,
        but does not check whether the valuation of the discriminant
        is divisible by $l$.

        If the elliptic curve has additive reduction, will raise a
        ValueError since the trace of Frobenius is not well-defined in
        that case.

        INPUT:

        - ``prime`` -- A prime of the definition field of this Frey
          curve. This must be a prime number if the definition field
          is $\QQ$ or a maximal ideal of the corresponding ring of
          integers otherwise.

        - ``power`` -- A strictly positive integer (default: 1). If
          set to a value higher than 1 will compute the trace of the
          Frobenius element to the given power instead of the
          Frobenius element itself.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        The algebraic number that is the Frobenius element at the
        given prime to the given power under the $l$-adic or mod $l$
        reprentation assuming that they are unramified. If this would
        depend on the parameters will return a conditional value of
        possible values instead.

        """
        pAdics = pAdicBase(self.definition_field(), prime)
        if power > 1:
            D = len(pAdics.residue_field())
            T = self.trace_of_frobenius(prime, condition=condition,
                                        verbose=verbose,
                                        precision_cap=precision_cap)
            result = self._power_trace_formula(power)
            return apply_to_conditional_value(lambda t: result(t, D), T)
        if condition is None:
            condition = self._condition
        red_type = self.reduction_type(prime, condition=condition,
                                       verbose=verbose,
                                       precision_cap=precision_cap)
        result = dict()
        if isinstance(red_type, ConditionalValue):
            for val, con in red_type:
                for a, T in self._trace_of_frobenius(pAdics, val, con, verbose,
                                                     precision_cap).items():
                    if a in result:
                        result[a] = result[a].merge(T)
                    else:
                        result[a] = T
        else:
            for a, T in self._trace_of_frobenius(pAdics, red_type, condition,
                                                 verbose,
                                                 precision_cap).items():
                if a in result:
                    result[a] = result[a].merge(T)
                else:
                    result[a] = T
        
        if len(result) == 1:
            return list(result)[0]
        else:
            Tls = [(a, pAdicTree(variables=self.parameters(), root=T))
                   for a, T in result.items()]
            return ConditionalValue([(a, TreeCondition(T)) for a, T in Tls])

    @cached_method
    def _power_trace_formula(self, n):
        r"""Give the formula to compute the trace of a matrix power.

        Given a 2-by-2 matrix $A$, the trace of $A^n$ for some $n \ge
        1$ can be expressed in terms of the trace and determinant of
        $A$ with a formula. This function gives this formula.

        """
        R.<x,y> = QQ[]
        f = x^n + y^n
        return polynomial_to_symmetric(f)

    @cached_method(key=lambda self, add, c, alg, v, prec, path:
                   ((self._condition if c is None else c),
                    (tuple(self.primes_of_possible_additive_reduction())
                     if add is None else tuple(add)), prec))
    def newform_candidates(self, bad_primes=None, condition=None,
                           algorithm='sage', verbose=False, precision_cap=20,
                           path=None):
        r"""Compute newforms that could be associated to this Frey curve.

        Given a Frey curve defined over the rationals, modularity
        tells us that its $l$-adic galois representations are
        isomorphic those arising from some modular form of level equal
        to the conductor of this curve.

        If for all prime numbers $p$ except for those in a finite set
        $S$ the order of $p$ in the discriminant is divisible by $l$
        and the conductor exponent at $p$ is at most $1$, level
        lowering results tell us that the same is true for the mod $l$
        representation, but that the associated newform in this case
        has level equal to the part of the conductor containing the
        primes in $S$. This function will assume we are in this case
        and return all the newforms of this level.

        INPUT:
        
        - ``bad_primes`` -- A list of prime numbers or None (default:
          None). This should be the list of prime numbers $p$ for
          which the order of $p$ in the discriminant is not divisible
          by $l$ or for which the curve has additive reduction. If set
          to None will be initialized as the result of
          :meth:`primes_of_possible_additive_reduction`, which might
          only contain the prime numbers at which the curve has
          additive reduction.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``algorithm`` -- The algorithm that should be used to
          compute the newforms. For possible options look at the
          function :func:`get_newforms`.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        - ``path`` -- An argument that might be required if the
          argument algorithm is set to a particular value. See the
          function :func:`get_newforms` for more explanation.

        OUTPUT:

        A list of newforms $f$ for which the mod $l$ galois
        representation of $f$ and this Frey curve might be
        isomorphic. If this list depends on the values of the
        parameters, returns a conditional value containing the
        possible lists with the associated conditions on the
        parameters.

        """
        if self.definition_field() != QQ:
            raise ValueError("Can only find newforms associated to " +
                             "Frey curves over the rationals.")
        if bad_primes is None and verbose >= 0:
            print("Warning: The bad primes chosen by default only take into "+
                   "account primes of additive reduction.")
        N = self.conductor(additive_primes=bad_primes, condition=condition,
                           verbose=verbose, precision_cap=precision_cap).left()
        if isinstance(N, ConditionalExpression):
            N = N.value()
        if condition is None:
            condition = self._condition
        return apply_to_conditional_value(lambda level:
                                          get_newforms(level,
                                                       algorithm=algorithm,
                                                       path=path), N)

    def newforms(self, condition=None, bad_primes=None, algorithm='sage',
                 primes=50, verbose=False, precision_cap_conductor=20,
                 precision_cap_reduction=1, path=None):
        r"""Use :meth:`newform_candidates` and :func:`eliminate_by_trace`
        instead.

        INPUT:

        - ``condition`` -- A Condition giving the restrictions on the
          parameters on this Frey curve that should be considered. By
          default this will be set to the condition associated to this
          FreyCurve.

        - ``bad_primes`` -- An iterable containing prime ideals or
          prime numbers, if the field of definition is QQ, that
          contains all the primes at which this curve can have
          additive reduction. If set to None will compute this by
          using the method primes_of_possible_additive_reduction

        - ``algorithm`` -- One of the following values 'sage' -- to
          use sage to compute newforms (default) 'magma' -- to use
          magma to compute newforms

        - ``primes`` -- A list of prime numbers or a strictly positive
          integer (default: 50). This list gives all the primes at
          which the traces of frobenius of the different galois
          representations should be compared. If set to strictly
          positive integer, will be initialized as the list of all
          prime numbers less than the given number.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings.  If this method calls any
          method that accepts an argument verbose will pass this
          argument to it.  If such a method fulfills a minor task
          within this method and the argument verbose was larger than
          0, will instead pass 1 less than the given argument. This
          makes it so a higher value will print more details about the
          computation than a lower one.

        - ``precision_cap_conductor`` -- A strictly positive integer
          (default: 20) giving the maximal precision level to be used
          in p-Adic arithmetic when computing the conductor.

        - ``precision_cap_reduction`` -- A strictly positive integer
          (default: 1) giving the maximal precision level to be used
          in the p-Adic arithmetic when computing the reduction type
          at a given prime. Since this will do a computation for every
          prime lower than prime_cap, this might get very
          computational intensive if set to a value larger than 1.

        - ``path`` -- A string or None (default: None). A parameter
          only used if the algorithm is set to file, in which case
          this should be a path to the file from which to load
          newforms.

        OUTPUT:

        A list consisting of pairs with as first entry a newform that
        has a mod-l representation that has traces of frobenius that
        could match the traces of frobenius of this curve for all
        given primes, different from l. The second entry is an integer
        divisible by all prime numbers l for which this can be true.
        
        If the level of the newform might depend on a choice of
        parameters will instead give a conditional value wherein each
        value is of the form above and each condition corresponds to a
        single possible level.

        """
        if condition is None:
            condition = self._condition
        if bad_primes is None:
            bad_primes = self.primes_of_possible_additive_reduction()
        newforms = self.newform_candidates(bad_primes=bad_primes,
                                           condition=condition,
                                           algorithm=algorithm,
                                           precision_cap=
                                           precision_cap_conductor,
                                           verbose=verbose, path=path)
        if primes in ZZ and primes > 0:
            primes = prime_range(primes)
        primes = list(primes)
        for P in bad_primes:
            if P in ZZ and P in primes:
                primes.remove(P)
            elif P.smallest_integer() in primes:
                primes.remove(P.smallest_integer())
        return eliminate_by_traces(self, newforms, condition=condition,
                                   primes=primes,
                                   precision_cap=precision_cap_reduction,
                                   verbose=(verbose - 1 if verbose > 0
                                            else verbose))
    def _repr_(self):
        """Give a string representation of a Frey curve.

        .. NOTE:

        This is a direct copy from the code included
        in :class:`EllipticCurve_number_field`

        """
        b = self.ainvs()
        a = [z._coeff_repr() for z in b]
        s = "Frey curve defined by "
        s += "y^2 "
        if a[0] == "-1":
            s += "- x*y "
        elif a[0] == '1':
            s += "+ x*y "
        elif b[0]:
            s += "+ %s*x*y "%a[0]
        if a[2] == "-1":
            s += "- y "
        elif a[2] == '1':
            s += "+ y "
        elif b[2]:
            s += "+ %s*y "%a[2]
        s += "= x^3 "
        if a[1] == "-1":
            s += "- x^2 "
        elif a[1] == '1':
            s += "+ x^2 "
        elif b[1]:
            s += "+ %s*x^2 "%a[1]
        if a[3] == "-1":
            s += "- x "
        elif a[3] == '1':
            s += "+ x "
        elif b[3]:
            s += "+ %s*x "%a[3]
        if a[4] == '-1':
            s += "- 1 "
        elif a[4] == '1':
            s += "+ 1 "
        elif b[4]:
            s += "+ %s "%a[4]
        s = s.replace("+ -","- ")
        s += "over %s "%(self.definition_field(),)
        s += "with parameters %s"%(self.parameters(),)
        return s

class FreyQcurve(FreyCurve, Qcurve):
    r"""A Frey-Hellegouarch curve that is also a Q-curve.

    .. SEE_ALSO::

        :class:`FreyCurve`
        :class:`Qcurve`

    """
    def __init__(self, curve, parameter_ring=ZZ, condition=None, **kwds):
        r"""Initializes a Frey Q-curve.

        This initialization calls the initialization of both
        :class:`Qcurve` and :class:`FreyCurve`. Note however that for
        the initialization of the first the parameter
        `guessed_degrees` is always set to the empty list as there is
        no way implemented to guess isogenies of a Frey curve of a
        given degree.

        INPUT:

        - ``curve`` -- An elliptic curve or any argument that would
          produce such a curve when passed to the constructor
          :func:`EllipticCurve`. The elliptic curve should be defined
          over some (multivariate) polynomial ring $R$ which in turn
          is defined over some Galois number field $L$.

        - ``parameter_ring`` -- The ring of integers of a subfield $K$
          of $L$ (default: ZZ). This is the ring in which the
          variables of the polynomial ring over which this curve is
          defined can take values.

        - ``condition`` -- An instance of :class:`Condition` or None
          (default: None) giving a condition which must hold for the
          values of the variables of $R$. If set to None will assume
          that all values for these variables are allowed.

         - ``isogenies`` -- A dictionary (default: {}) with as keys
           elements s of the Galois group of the base field of the
           Q-curve and as values the corresponding isogeny from this
           curve to the Galois conjugate by s. Such isogenies must be
           defined over the ring $R$ over which the curve is defined
           or a polynomial ring in the same variable defined over an
           extension of the field $L$. The isogeny can be given as
           either an isogeny as a Sage object; a tuple of a rational
           function in $x$ and an algebraic number, that are
           respectively the $x$-coordinate map of the isogeny and the
           induced scalar multiplication on the differentials; or a
           tuple of three rational functions $F$, $G$, $H$ in $x$ such
           that the isogeny is $(x, y) \mapsto (F(x), G(x) y + H(x))$
           outside points mapping to infinity.

        """
        FreyCurve.__init__(self, curve, parameter_ring=parameter_ring,
                           condition=condition)
        Qcurve.__init__(self, curve, **kwds)

    def _init_curve(self, curve):
        r"""Initialize the underlying elliptic curve.

        This overwrites the method found in :class:`Qcurve`, such that
        the methods of that class can work with this curve as if it
        was defined over its definition field, rather than the
        polynomial ring that is its base ring.

        When this method is called most things have already been
        initialized by the initialization from :class:`FreyCurve`, but
        we initialize again to make sure the definition field will be
        galois over $\QQ$.

        """
        K = self.definition_field()
        if not is_NumberField(K) or not K.is_galois():
            raise ValueError("The ring " + str(K) +
                             " is not a Galois number field.")

    @cached_method(key=lambda self, sigma, change : (self._galois_cache_key(sigma), change))
    def galois_conjugate(self, sigma, change_ring=None):
        r"""Give the Galois conjugate of this curve.

        INPUT:

        - ``sigma`` -- A Galois homomorphism of some number field

        - ``change_ring`` -- A field homomorphism from the definition
          field of this curve to another field or None (default). If
          set to a value other than None, this function will return
          the Galois conjugate over the field that is the codomain of
          this homomorphism.

        OUTPUT:
        
        The galois conjugate of this curve by the galois homomorphism
        which extends to a common galois homomorphism over the
        algebraic closure of Q as sigma. This will be an elliptic
        curve and not a Q-curve

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?
            sage: sigma = K.galois_group().gens()[0]
            sage: E.galois_conjugate(sigma)
            Elliptic Curve defined by y^2 = x^3 + 12*x^2 + (-18*t+18)*x over Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?

        """
        sigma = galois_field_change(sigma, self.definition_field())
        change = sigma.as_hom()
        if change_ring != None:
            change = _concat_maps(change, change_ring)
        R = self.base_ring()
        S = R.change_ring(change.codomain())
        change = R.Hom(S)(change)
        return EllipticCurve(self.a_invariants()).change_ring(change)

    def definition_field(self):
        r"""Give the field over which this Frey curve is defined.

        Even though the Frey curve is defined over some (multivariate)
        polynomial ring $R$ over some number field $L$, since the
        variables of $R$ are assumed to have values in some subfield
        $K$ of $L$ the curve can be assumed to be defined over $L$.

        OUTPUT:
        
        The base ring of the polynomial ring over which this Frey
        curve is defined.

        .. SEE_ALSO::

            :meth:`base_ring`,
            :meth:`parameters`

        """
        return self.base_ring().base()

    def _get_isogeny_field(self, sigma):
        r"""Get the field over which an isogeny is defined"""
        return self._phi_x[sigma].base_ring().base_ring()

    def _get_isogeny(self, sigma, change=None):
        r"""Give the x and y rational maps of an isogeny

        INPUT:

        - ``sigma`` -- A galois homomorphism of the definition field
          of this curve.

        - ``change`` -- A field homomorphism or None (default). If not
          None, should be a field homomorphism from the field over
          which the isogeny is currently defined.

        OUTPUT:

        A tuple consisting of the x-map and the y-map of the isogeny
        from this curve to the `sigma` conjugate. If change was set to
        a field homomorphism, these maps will be defined over the
        codomain of this map.

        """
        phi_x = self._phi_x[sigma]
        phi_y = self._phi_y[sigma]
        if change != None:
            change = _write_as_im_gen_map(change)
            Rxy = phi_y.parent().base()
            R = Rxy.base_ring()
            S = R.change_ring(change.codomain())
            Sxy = Rxy.change_ring(S).fraction_field()
            change = R.Hom(S)(change)
            phi_x = (Sxy(phi_x.numerator().change_ring(change)) /
                     Sxy(phi_x.denominator().change_ring(change)))
            phi_y = (Sxy(phi_y.numerator().change_ring(change)) /
                     Sxy(phi_y.denominator().change_ring(change)))
        return phi_x, phi_y

    def _add_isogenies_of_degree(self, degree, verbose=False):
        r"""Attempt to find isogenies of a given degree"""
        K = self.definition_field()
        R = self.base_ring()
        F = R.fraction_field()
        G = K.galois_group()
        S.<x> = K[]
        fd = self.torsion_polynomial(degree)
        E = self.galois_conjugate(G.identity()).change_ring(F)
        sF = {s : F.Hom(F)(R.Hom(R)(s.as_hom())) for s in G}
        Es = {s : E.change_ring(sF[s]) for s in G}
        js = {s : Es[s].j_invariant() for s in G}
        for g, e in fd.factor():
            if g.degree() != 1:
                continue
            g = g.change_ring(F) / F(g[1])
            psi = E.isogeny(g)
            Ed = psi.codomain()
            jd = Ed.j_invariant()
            for sigma in G:
                if js[sigma] != jd:
                    continue
                c4s, c6s = Es[sigma].c_invariants()
                c4d, c6d = Ed.c_invariants()
                if jd == 0:
                    m, um = 6, c6d/c6s
                elif jd == 1728:
                    m, um = 4, c4d/c4s
                else:
                    m, um = 2, (c6d*c4s)/(c6s*c4d)
                if um.numerator().monomials() != um.denominator().monomials():
                    # Can not be an element of K
                    continue
                m0 = um.numerator().monomials()[0]
                um0 = (um.numerator().monomial_coefficient(m0) /
                       um.denominator().monomial_coefficient(m0))
                if um.numerator() != um0 * um.denominator():
                    # Not an element of K
                    continue
                fu = x^m - um0
                for gu, e in fu.factor():
                    L.<lu> = K.extension(gu)
                    RL = R.change_ring(L)
                    FL = RL.fraction_field()
                    iota = F.Hom(FL)(R.Hom(RL)(_write_as_im_gen_map(K.hom(L))))
                    EL = E.change_ring(iota)
                    EdL = Ed.change_ring(iota)
                    EsL = Es[sigma].change_ring(iota)
                    s = (lu*EsL.a1() - EdL.a1())/2
                    r = (lu^2*EsL.a2() - EdL.a2() + s^2 + s*EdL.a1())/3
                    t = (lu^3*EsL.a3() - EdL.a3() + r*EdL.a1())
                    if (lu^4*EsL.a4() != (EdL.a4() - s*EdL.a3() + 2*r*EdL.a2()
                                           - (t + r*s)*EdL.a1() + 3*r^2 - 2*s*t) or
                        lu^6*EsL.a6() != (EdL.a6() + r*EdL.a4() + r^2*EdL.a2() +
                                           r^3 - t*EdL.a3() - t^2 - r*t*EdL.a1())):
                        continue
                    psi = EL.isogeny(g.change_ring(iota))
                    phi = WeierstrassIsomorphism(E=EdL, urst=(lu, r, s, t))
                    psi.set_post_isomorphism(phi)
                    if verbose > 0:
                        print("Degree %s isogeny found for"%degree, sigma)
                    self._add_isogeny(sigma, psi)
    
    @cached_method(key=lambda self, sigma, change :
                   (self._galois_cache_key(sigma), change))
    def isogeny_x_map(self, sigma, change_ring=None):
        r"""Return the x-coordinate rational map of the isogeny from this curve
        to a Galois conjugate.

        The x-coordinate of the image of a point under an isogeny can
        be described as a rational function of the x-coordinate of the
        corresponding point in the domain.

        INPUT:
        
        - ``sigma`` -- A Galois homomorphism of a number field

        - ``change_ring`` -- A field homomorphism from the definition
          field of this Q-curve or None (default). If set to a value
          other than None, the base field of the returned rational map
          will be changed by this homomorphism.

        OUTPUT:

        A rational function in $x$ over the definition field of this
        Q-curve that gives the $x$-coordinate of an image point of the
        isogeny as a rational function in the $x$ coordinate of the
        origin point. The isogeny is the registered isogeny from this
        curve to the `sigma` Galois conjugate of this curve. If
        `change_ring` was set to None, the returned rational function
        will be defined over the codomain of `change_ring`.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-1)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: G.<s> = K.galois_group()
            sage: E.isogeny_x_map(s)
            (-1/2*x^2 - 6*x - 9*t - 9)/x
            sage: E.isogeny_x_map(s^2)
            x

        TESTS::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-1)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: G.<s> = K.galois_group()
            sage: F = E.isogeny_x_map(s)
            sage: F.parent()
            Fraction Field of Univariate Polynomial Ring in x over Number Field in t with defining polynomial x^2 + 1 with t = 1*I
            sage: F.parent().base_ring() == E.definition_field()
            True

        .. SEEALSO::

            :meth:`isogeny_scalar`,
            :meth:`definition_field`

        """
        sigma = galois_field_change(sigma, self.definition_field())
        R = self.base_ring()
        Rx.<x> = PolynomialRing(R)
        Fx = Rx.fraction_field()
        phi_x = self._phi_x[sigma]
        num = phi_x.numerator()
        S = num.parent()
        den = phi_x.denominator()
        _, iota = write_as_extension(self._to_Kphi, give_map=True)
        result = (Fx(sum(sum(iota(cf).list()[0] * R(m)
                             for cf, m in num.monomial_coefficient(S.monomial(*e)))
                         * x^e[0] for e in num.exponents())) /
                  Fx(sum(sum(iota(cf).list()[0] * R(m)
                             for cf, m in den.monomial_coefficient(S.monomial(*e)))
                         * x^e[0] for e in den.exponents())))
        if change_ring is None:
            return result
        R2 = R.change_ring(change_ring.codomain())
        Rx2 = Rx.change_ring(R2)
        Fx2 = Rx2.fraction_field()
        change_ring = Fx.Hom(Fx2)(Rx.Hom(Rx2)(R.Hom(R2)(change_ring)))
        return change_ring(result)

    def _isogeny_data(self, iota):
        r"""Give the isogeny data of this curve over a given field.
            
        INPUT:

        - ``iota`` -- A field homomorphism from the definition field
          of this Q-curve to another Galois number field $K$. 

        OUTPUT:

        A dictionary of which the keys are the elements of the Galois
        group of $K$, and the value for each element sigma is a tuple
        of a rational function in $x$ and an algebraic number, that
        are respectively the $x$-coordinate map of the isogeny from
        the sigma conjugate of this curve to this curve and the
        induced scalar multiplication on the differentials.

        """
        G = iota.codomain().galois_group()
        F = self.isogeny_x_map
        l = self.isogeny_scalar
        Kphi = self.complete_definition_field()
        _, to_L, from_Kphi = composite_field(iota, self._to_Kphi,
                                             give_maps=True)
        iota = _concat_maps(iota, to_L)
        L, to_L = write_as_extension(to_L, give_map=True)
        from_Kphi = _concat_maps(from_Kphi, to_L)
        iota = _concat_maps(iota, to_L)
        R = self.base_ring().change_ring(L)
        iota = self.base_ring().Hom(R)(iota)
        Rx = PolynomialRing(R, names=["x"]).fraction_field()
        x = Rx.gens()[0]
        return {s : (Rx(F(s).numerator().change_ring(iota) /
                        F(s).denominator().change_ring(iota)),
                     l(s).change_ring(from_Kphi))
                for s in G}

    def base_extend(self, R):
        result = FreyCurve.base_extend(self, R)
        if (isinstance(result, FreyCurve) and
            is_NumberField(result.definition_field())):
            K = self.definition_field()
            L = result.definition_field()
            if (hasattr(R, "codomain") and
                R.domain() == K and
                R.codomain() == L):
                K_to_L = R
            else:
                K_to_L = K.embeddings(L)
                if len(K_to_L) > 0:
                    K_to_L = K_to_L[0]
                else:
                    K_to_L = None
            if K_to_L != None:
                return FreyQcurve(result,
                                  isogenies=self._isogeny_data(K_to_L),
                                  parameter_ring=result._R,
                                  condition=result._condition)
        return result

    def _minimal_a_invariants(self, from_min):
        r"""Give the a_invariants over the minimal definition_field

        INPUT:

        - ``from_min`` -- The map from the minimal definition field to
          the definition field.

        """
        _, iota = write_as_extension(from_min, give_map=True)
        R = self.base_ring().change_ring(from_min.domain())
        return [sum(iota(a.monomial_coefficient(m)).list()[0] * R(m)
                    for m in a.monomials()) for a in self.a_invariants()]
        
    def _minimal_isogeny_x_maps(self, Kmin, from_min, min_map):
        r"""Give the x-coordinate maps of isogenies for :meth:minimize_fields

        INPUT:

        - ``Kmin`` -- The minimal definition field
        
        - ``from_min`` -- The map from the minimal definition field to
          the definition field.

        - ``min_map`` -- The map from the minimal definition field to
          the minimal complete definition field.

        """
        result = {}
        G = Kmin.galois_group()
        _, iota = write_as_extension(from_min, give_map=True)
        for s in G:
            Fs = self.isogeny_x_map(s)
            Fsnum = Fs.numerator()
            Fsden = Fs.denominator()
            R = Fsnum.parent().base().change_ring(Kmin)
            S = Fsnum.parent().change_ring(R)
            Fsnum = sum(sum(iota(cf.monomial_coefficient(mR)).list()[0] *
                            R(mR) for mR in cf.monomials()) *
                            S(mS) for cf, mS in
                        [(Fsnum.monomial_coefficient(mS), mS)
                         for mS in Fsnum.monomials()])
            Fsden = sum(sum(iota(cf.monomial_coefficient(mR)).list()[0] *
                            R(mR) for mR in cf.monomials()) *
                            S(mS) for cf, mS in
                        [(Fsden.monomial_coefficient(mS), mS)
                         for mS in Fsden.monomials()])
            R2 = R.change_ring(min_map.codomain())
            S2 = S.change_ring(R2)
            change = S.Hom(S2)(R.Hom(R2)(min_map))
            Fsnum = change(Fsnum)
            Fsden = change(Fsden)
            result[s] = S2.fraction_field()(Fsnum / Fsden)
        return result

    def minimize_fields(self, names=None):
        r"""Attempt to minimize the fields associated to this curve.

        INPUT:

        - ``names`` -- A tuple or list of two strings or None
          (default). These two strings will be respectively the name
          of the variable of the definition field and the complete
          definition field of the returned Q-curve. If set to None
          will be initialized as the names of the definition field and
          complete definition field of this curve.

        OUTPUT:

        A Q-curve isomorphic to this Q-curve, but with a definition
        field, complete definition field and decomposition field that
        are as small as possible. The isomorphism is given by
        compatible embeddings of the mentioned fields.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: iota = K.embeddings(CyclotomicField(12))[0]
            sage: E2 = E.change_ring(iota)
            sage: E3 = E2.minimize_fields(names=["t", "s"])
            sage: E.definition_field()
            Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?
            sage: E2.definition_field()
            Cyclotomic Field of order 12 and degree 4
            sage: E3.definition_field()
            Number Field in zeta120 with defining polynomial x^2 - 3 with zeta120 = 1.732050807568878?
            sage: E.complete_definition_field()
            Number Field in lu with defining polynomial x^4 - 2*x^2 + 25
            sage: E2.complete_definition_field()
            Number Field in zeta12lu with defining polynomial x^8 - 18*x^6 + 239*x^4 - 1638*x^2 + 6241
            sage: E3.complete_definition_field()
            Number Field in s with defining polynomial x^4 - 38*x^2 + 1225

        """
        Kmin, Kphi, from_min, min_map, phi_map = self._minimize_fields(names=names)
        ainvs = self._minimal_a_invariants(from_min)
        F = self._minimal_isogeny_x_maps(Kmin, from_min, min_map)
        l = self._minimal_isogeny_scalars(Kmin, phi_map)
        isogenies = {s : (F[s], l[s]) for s in Kmin.galois_group()}
        return FreyQcurve(ainvs, isogenies=isogenies,
                          parameter_ring=self._R,
                          condition=self._condition)

    def twist(self, gamma):
        r"""Give the twist of this Frey Q-curve by a given element gamma.

        If this curve was given by .. MATH::

            E : y^2 = x^3 + a_2 x^2 + a_4 x + a_6

        the twisted curve is given by .. MATH::
        
            E : y^2 = x^3 + \gamma a_2 x^2 + \gamma^2 a_4 x
                      + \gamma^3 a_6

        INPUT:

        - ``gamma`` -- An element of a number field.

        OUTPUT:
        
        A Frey Q-curve which is the twist of this Q-curve by
        gamma. The definition field of this new curve will be the
        smallest possible field over which it is completely defined as
        a Q-curve.

        """
        ainvs, isogenies = self._twist(gamma)
        return FreyQcurve(ainvs, isogenies=isogenies,
                          parameter_ring=self._R,
                          condition=self._condition).minimize_fields()
    
    def conductor_restriction_of_scalars(self, additive_primes=None,
                                         condition=None, verbose=False,
                                         precision_cap=20):
        r"""Give the conductor of the restriction of scalars of this Frey
        Q-curve.

        INPUT:

        - ``additive_primes`` -- A list containing primes of the
          definition field of this curve or None (default: None). The
          primes in this list should be given as prime number if the
          definition field is $\QQ$ or as maximal ideals
          otherwise. This list should include all the primes at which
          this curve could have additive reduction. If set to None
          will compute this by using the method
          :meth:`primes_of_possible_additive_reduction`.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey Q-curve
          should satisfy. If set to None will use the condition stored
          in this Frey Q-curve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings.  A higher value will cause
          more messages to be printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        The conductor of the restriction of scalars of this curve over
        the decomposition field. This will be a conditional expression
        containing on the left side a (conditional) expression of the
        part of the conductor coming from primes in `additive_primes`,
        whilst the right hand side is a string describing how to
        compute the part of the conductor coming from primes coprime
        to the primes in `additive_primes`. The latter contains the
        operator Rad_P which refers to taking the radical of an
        expression ignoring those primes in additive_primes.

        """
        K0 = self.definition_field()
        K = self.decomposition_field()
        if K0 != K:
            iota = _concat_maps(self._to_Kphi, self._to_Kdec)
            E = self.change_ring(iota)
        else:
            E = self
        if additive_primes is None:
            additive_primes = copy(E.primes_of_possible_additive_reduction())
            for p in K.discriminant().prime_factors():
                for P in K.primes_above(p):
                    if P.ramification_index() > 1 and P not in additive_primes:
                        additive_primes.append(P)
        # Proposition 1 of Milne, On the arithmetic of Abelian varieties
        N = E.conductor(additive_primes=additive_primes,
                        condition=condition,
                        verbose=verbose,
                        precision_cap=precision_cap)
        additive_part = N.left()
        Dsqr = K.discriminant()^2
        if isinstance(additive_part, ConditionalExpression):
            additive_factors = N.left().factors()
            left_factors = {}
            for f in additive_factors:
                for p,e in f.absolute_norm().factor():
                    if e == 1:
                        e = additive_factors[f]
                    else:
                        e = e * additive_factors[f] 
                    if p in left_factors and left_factors[p] != 0:
                        left_factors[p] = left_factors[p] + e
                    elif e != 0:
                        left_factors[p] = e
            disc_factors = Dsqr.factor()
            for p, e in disc_factors:
                if p in left_factors and left_factors[p] != 0:
                    left_factors[p] = left_factors[p] + e
                elif e != 0:
                    left_factors[p] = e
            if hasattr(disc_factors, 'unit') and disc_factors.unit() != 1:
                left = (Dsqr.factor().unit() *
                        product(p^e for p, e in left_factors.items()))
            else:
                left = product(p^e for p,e in left_factors.items())
        else:
            left = additive_part.absolute_norm() * Dsqr
        return ConditionalExpression(N.operator(),
                                     left,
                                     "Norm(" + N.right() +")")

    def newform_levels(self, bad_primes=None, condition=None, verbose=False,
                       precision_cap=20):
        r"""Compute the levels of newforms that could be associated to this
        Frey Q-curve.

        Each non-CM Q-curve is the quotient of a $\Q$-simple variety
        of GL_2-type, which in turn is isogenous to an abelian
        varietyr associated to a newform. The $\lambda$-adic galois
        representation of this newform is isomorphic to the $l$-adic
        galois representation of the Q-curve when restricted to a
        common subgroup of the absolute galois group of $\QQ$. Here
        $\lambda$ is a prime dividing $l$ in the coefficient field of
        the newform.

        The conductor of an abelian variety associated to a newform is
        $N^n$, where $N$ is the level of the newform and $n$ is the
        dimension of the variety. If the Q-curve decomposes, the
        factors of its restriction of scalars form abelian varieties
        of associated newforms. These newforms are directly related to
        the splitting maps of the Q-curves, in the sense that they are
        twists of one another by the inverse of the twist characters
        and their characters are the inverse of the splitting
        characters. Using results about the change in level when
        twisting a newform and the conductor of the restriction of
        scalars, a guess for the levels of the newforms can be made.

        Let $E$ be an elliptic and $P$ be a prime of its decomposition
        field for which the order of $P$ in the discriminant of $E$ is
        divisible by a prime number $l$ and for which $E$ does not
        have additive reduction. In that case is the mod $l$ galois
        representation of $E$ unramified at $P$. In case $E$ is a
        $\Q$-curve and $P$ is not ramified in the decomposition field,
        this implies that the corresponding mod $\lambda$
        representation of an associated newform is unramified at the
        prime number $p$ below $P$. In this case we would be able to
        find a newform of a lower level that has an isomorphic mod
        $\lambda$ representation. To be precise the lower level is the
        part of the level that is coprime to $p$.

        In the case we have a Frey Q-curve $E$ and only a finite set
        $S$ of bad primes of its decompisition field, i.e. primes $P$
        for which the curve $E$ has additive reduction, that ramify in
        the decomposition field or for which their order in the
        discriminat of $E$ is not divisible by $l$, we know by level
        lowering that the mod $l$ representation of $E$ is isomorphic
        to the mod $\lambda$ galois representation of newforms with a
        level only divisible by prime numbers below primes in the set
        $S$. This function computes the possible levels for a given
        set $S$ of bad primes.

        INPUT:
        
        - ``bad_primes`` -- A list of primes of the decomposition
          field of this curve or None (default: None). These primes
          should be given as prime numbers if the decomposition field
          is $\QQ$ or as prime ideal otherwise. This should be the
          list of all the bad primes, i.e. primes for which this curve
          has additive reduction, primes that ramify in the
          decomposition field and primes of which the order in the
          disriminant of this curve is not divisible by $l$. If set to
          None will be initialized as the result of
          :meth:`primes_of_possible_additive_reduction` together with
          all primes that ramify in the decomposition field, which
          might omit some primes for which the discriminant of this
          curve is not an $l$-th power.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        A list of tuples, each tuple representing one of the options
        for the levels of the newforms associated to this Q-curve. The
        $i$-th entry of such a tuple is the level of a newform
        corresponding to the $i$-th conjugacy class of splitting maps,
        as returned by :meth:`splitting_map`.

        If the outcome depends on the values of the parameters, will
        return a conditional value containing such values with the
        appropiate conditions on the parameters.

        """
        if condition is None:
            condition = self._condition
        N = self.conductor_restriction_of_scalars(additive_primes=bad_primes,
                                                  condition=condition,
                                                  verbose=verbose,
                                                  precision_cap=precision_cap)
        N = N.left()
        if isinstance(N, ConditionalExpression):
            N = N.value()
        return apply_to_conditional_value(lambda Ni:
                                          Qcurve.newform_levels(self, N=Ni), N)

    @cached_method(key=lambda self, add, c, alg, prec, v, path:
                   ((self._condition if c is None else c),
                    (tuple(self.primes_of_possible_additive_reduction())
                     if add is None else tuple(add)),
                    prec))
    def newform_candidates(self, bad_primes=None, condition=None,
                           algorithm='sage', verbose=False, precision_cap=20,
                           path=None):
        r"""Compute newforms that could be associated to this Frey Q-curve.

        Each non-CM Q-curve is the quotient of a $\Q$-simple variety
        of GL_2-type, which in turn is isogenous to an abelian
        varietyr associated to a newform. The $\lambda$-adic galois
        representation of this newform is isomorphic to the $l$-adic
        galois representation of the Q-curve when restricted to a
        common subgroup of the absolute galois group of $\QQ$. Here
        $\lambda$ is a prime dividing $l$ in the coefficient field of
        the newform.

        The conductor of an abelian variety associated to a newform is
        $N^n$, where $N$ is the level of the newform and $n$ is the
        dimension of the variety. If the Q-curve decomposes, the
        factors of its restriction of scalars form abelian varieties
        of associated newforms. These newforms are directly related to
        the splitting maps of the Q-curves, in the sense that they are
        twists of one another by the inverse of the twist characters
        and their characters are the inverse of the splitting
        characters. Using results about the change in level when
        twisting a newform and the conductor of the restriction of
        scalars, a guess for the levels of the newforms can be made.

        Let $E$ be an elliptic and $P$ be a prime of its decomposition
        field for which the order of $P$ in the discriminant of $E$ is
        divisible by a prime number $l$ and for which $E$ does not
        have additive reduction. In that case is the mod $l$ galois
        representation of $E$ unramified at $P$. In case $E$ is a
        $\Q$-curve and $P$ is not ramified in the decomposition field,
        this implies that the corresponding mod $\lambda$
        representation of an associated newform is unramified at the
        prime number $p$ below $P$. In this case we would be able to
        find a newform of a lower level that has an isomorphic mod
        $\lambda$ representation. To be precise the lower level is the
        part of the level that is coprime to $p$.

        In the case we have a Frey Q-curve $E$ and only a finite set
        $S$ of bad primes of its decompisition field, i.e. primes $P$
        for which the curve $E$ has additive reduction, that ramify in
        the decomposition field or for which their order in the
        discriminat of $E$ is not divisible by $l$, we know by level
        lowering that the mod $l$ representation of $E$ is isomorphic
        to the mod $\lambda$ galois representation of newforms with a
        level only divisible by prime numbers below primes in the set
        $S$. This function computes the possible levels and all
        newforms at those levels for this curve given such a set $S$.

        INPUT:
        
        - ``bad_primes`` -- A list of primes of the decomposition
          field of this curve or None (default: None). These primes
          should be given as prime numbers if the decomposition field
          is $\QQ$ or as prime ideal otherwise. This should be the
          list of all the bad primes, i.e. primes for which this curve
          has additive reduction, primes that ramify in the
          decomposition field and primes of which the order in the
          disriminant of this curve is not divisible by $l$. If set to
          None will be initialized as the result of
          :meth:`primes_of_possible_additive_reduction` together with
          all primes that ramify in the decomposition field, which
          might omit some primes for which the discriminant of this
          curve is not an $l$-th power.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``algorithm`` -- The algorithm that should be used to
          compute the newforms. For possible options look at the
          function :func:`get_newforms`.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        - ``path`` -- An argument that might be required if the
          argument algorithm is set to a particular value. See the
          function :func:`get_newforms` for more explanation.

        OUTPUT:

        A list of newforms $f$ for which the mod $l$ galosi
        representations of $f$ and this Frey curve might be
        isomorphic. If this list depends on the values of the
        parameters, returns a conditional value containing the
        possible lists with the associated conditions on the
        parameters.

        """
        if condition is None:
            condition = self._condition
        levels = self.newform_levels(bad_primes=bad_primes,
                                      condition=condition,
                                      verbose=(verbose - 1 if verbose > 0
                                               else verbose),
                                      precision_cap=precision_cap)
        return apply_to_conditional_value(lambda levelsi, con:
                                          self._newform_candidates(levelsi,
                                                                   con &
                                                                   condition,
                                                                   algorithm,
                                                                   path,
                                                                   verbose),
                                          levels, use_condition=True,
                                          default_condition=condition)

    def _newform_candidates(self, levels, condition, algorithm, path, verbose):
        r"""Implementation of :meth:`newform_candidates`

        INPUT:

        - ``levels`` -- List of possible levels for the newforms

        - ``condition`` -- The condition on the parameters to get
          these levels

        - ``algorithm`` -- The argument `algorithm`

        - ``path`` -- The argument ``path``

        - ``verbose`` -- Verbosity argument

        """
        result = []
        done_levels = []
        characters = [(eps^(-1)).primitive_character()
                      for eps in self.splitting_character('conjugacy')]
        KE = self.splitting_image_field()
        for levelsi in levels:
            level, eps = min(zip(levelsi, characters),
                             key=lambda x: x[0])
            if (level, eps) in done_levels:
                continue # Already computed, continue on with the next
            if verbose > 0:
                print("Computing newforms of level %s and character %s"%(level, eps))
            for f in get_newforms(level, character=eps,
                                  algorithm=algorithm, path=path):
                Kf = f.coefficient_field()
                K = composite_field(KE, Kf)
                K = K.galois_closure(names=K._names)
                for phi in Kf.embeddings(K):
                    f_phi = f.copy()
                    f_phi.set_embedding(K, phi)
                    result.append(f_phi)
            done_levels.append((level, eps))
        return result

    def trace_of_frobenius(self, prime, power=1, splitting_map=0,
                           condition=None, verbose=False, precision_cap=20):
        r"""Compute the trace of a Frobenius element under the Galois
        representation associated to this curve.

        This function computes an algebraic number that, when cast to
        the appropriate field, is equal to the trace of (a power of) a
        Frobenius element of the given prime number under the Galois
        representation associated to the given splitting map. It can
        only compute this number in case

        - This curve has good reduction at the prime number and the
          reduction of the isogeny associated to the Frobenius element
          is separable

        - This curve has multiplicative reduction at the prime number
          and the isogeny associated to the Frobenius element has
          square free degree.

        Given that this Q-curve decomposes over its decomposition
        field, one can associate to each splitting map $\beta$ a
        $\QQ$-simple abelian variety $A_\beta$ of $GL_2$-type which
        has this Q-curve as a 1-dimensional quotient. For each finite
        prime $\lambda$ in the image field of $\beta$, this $A_\beta$
        defines a 2-dimensional $\lambda$-adic Galois representation
        of the absolute Galois group of $\QQ$ that extends the
        $l$-adic Galois representation of this curve over the
        decomposition field, where $l$ is the prime number below
        $\lambda$.

        .. NOTE::

        This method does not check whether the associated Galois
        representation is unramified. If it is not, but the conditions
        mentioned before are satisfied, this method will produce an
        answer corresponding to a Frobenius element. Note that the
        trace of a Frobenius element is in that case not unique.

        Since the Galois representations over the rationals extend the
        Galois representations of this curve over the decomposition
        field, each Galois representation is unramified at a prime p
        if its corresponding restriction to the decomposition field is
        unramified and p is unramified in the decomposition field.

        INPUT:

        - ``prime`` -- A prime number for which the curve either has
          good reduction and the reduction of the isogeny associated
          to a frobenius element of this prime is separable, or for
          which the curve has multiplicative reduction.

        - ``power`` -- A strictly positive integer (default: 1). If
          set to a value higher than 1 will compute the trace of the
          Frobenius element to the given power instead of the
          Frobenius element itself.

        - ``splitting_map`` -- A non-negative integer smaller than the
          number of splitting maps associated to this Q-curve
          (default: 0). This indicates the splitting map associated to
          the Galois representation for which the trace of Frobenius
          should be computed.

        - ``condition`` -- A Condition or None (default: None) giving
          the condition that the parameters of this Frey curve should
          satisfy. If set to None will use the condition stored in
          this FreyCurve instead.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such comments.
          If set to any negative value will also prevent the printing of
          any warnings.  A higher value will cause more messages to be
          printed.

        - ``precision_cap`` -- A strictly positive integer (default:
          20) giving the maximal precision level to be used in p-Adic
          arithmetic for the parameters.

        OUTPUT:

        An algebraic number that, when cast to the correct field, is
        the trace of a Frobenius element at `prime` to the power
        `power` under the $\lambda$-adic or mod $\lambda$ Galois
        representation associated to the splitting map given by
        `splitting_map` for $\lambda$ not dividing `prime`. If this
        would depend on the parameters will return a conditional value
        of possible values instead.

        """
        if power > 1:
            T = self.trace_of_frobenius(prime,
                                        splitting_map=splitting_map,
                                        condition=condition,
                                        verbose=verbose,
                                        precision_cap=precision_cap)
            D = self.determinant_of_frobenius(prime,
                                              splitting_map=splitting_map)
            result = self._power_trace_formula(power)
            return apply_to_conditional_value(lambda t: result(t, D), T)
        
        if condition is None:
            condition = self._condition
        if not self.does_decompose():
            raise ValueError("This Q-curve must decompose for " +
                             "the Galois representations over the " +
                             "rationals to extend the Galois " +
                             "representations over the decomposition " +
                             "field.")
        K = self.decomposition_field()
        if self.decomposition_field() != K:
            return self._over_Kdec().trace_of_frobenius(prime, power=power,
                                                        splitting_map=splitting_map,
                                                        condition=condition,
                                                        verbose=verbose,
                                                        precision_cap=precision_cap)
        if prime.divides(K.absolute_discriminant()):
            print("Warning: The decomposition field is ramified " +
                  "at " + str(prime) + " and the Galois " +
                  "representations over the rationals are therefore " +
                  "probably not unramified at " + str(prime))
        beta = self.splitting_map(splitting_map)

        def compute_trace1(P, con1):
            def compute_trace2(red_type, con2):
                if red_type == None: # Good reduction
                    return self._trace_of_frob_good(P, beta,
                                                    con2,
                                                    verbose,
                                                    precision_cap)
                elif red_type == 1 or red_type == -1: # Multiplicative reduction
                    return self._trace_of_frob_mult(P, beta,
                                                    con2,
                                                    verbose,
                                                    precision_cap)
                else: # Additive reduction
                    return None
                
            red_type = self.reduction_type(P, condition=con1,
                                           verbose=verbose,
                                           precision_cap=precision_cap)
            return apply_to_conditional_value(compute_trace2, red_type,
                                              use_condition=True,
                                              default_condition=con1)

        result = conditional_over_values(compute_trace1,
                                         K.primes_above(prime),
                                         start_condition=condition)
        if isinstance(result, ConditionalValue):
            for val, con in result:
                if val == None:
                    raise ValueError("Trace of frobenius can not be " +
                                     "computed in case " + str(con) +
                                     " holds.")
        elif result == None:
            raise ValueError("Trace of frobenius can not be " +
                             "computed in case " + str(con) +
                             " holds.")
        return result

    def _trace_of_frob_good(self, P, beta, condition, verbose,
                            precision_cap):
        r"""Implementation of meth:`trace_of_frobenius` in case of good reduction"""
        # Setting up some stuff needed everywhere
        K = self.definition_field()
        G = K.galois_group()
        Frob = G.artin_symbol(P)
        sE = FreyCurve(self.galois_conjugate(Frob),
                       condition=self._condition)
        sP = Frob(P)
        psi = self._phi_x[Frob], self._phi_y[Frob]
        R = P.residue_field()
        Rx.<x> = R[]
        Rx = Rx.fraction_field()
        pAdicsK = pAdicBase(K, P)
        pAdics = pAdicsK.pAdics_below(self._R)
        ext = pAdicsK.extension_multiplicity(pAdics)

        # Setting up the isogeny between the right minimal
        # models. Note that compute_trace1 and compute_trace2 will be
        # called near the end once for each different value of Emin
        # and sEmin respectively with the corresponding
        # condition. They can be considered as loops here, but note
        # that the final return value will be a combination of all
        # return values in each iteration of these loops..
        def compute_trace1(min_data, con1):
            Emin = min_data[0]
            phi_min = _rational_maps_of_urst(*min_data[1][1])
            psi1 = psi[0](phi_min), psi[1](phi_min)
            def compute_trace2(smin_data, con2):
                sEmin = smin_data[0]
                sphi_min = _rational_maps_of_urst(*smin_data[1][0])
                phi = sphi_min[0](psi1), sphi_min[1](psi1)

                # Making a loop over the possible reductions of Emin
                results = {}
                T = con2.pAdic_tree(pAdics=pAdics).root()
                E_val = ZZ(min(pAdicsK.valuation(cf)
                               for a in Emin.a_invariants()
                               for cf in a.coefficients())) // ext
                sE_val = ZZ(min(pAdicsK.valuation(cf)
                                for a in sEmin.a_invariants()
                                for cf in a.coefficients())) // ext
                T_val = T.minimum_full_level()
                val = max(1, min(1 - E_val, 1 - sE_val, T_val,
                                 precision_cap))
                for node in T.children_at_level(val):
                    fill = Emin.base_ring().hom([K(xi) for xi in
                                                 node.representative()])
                    Ered = EllipticCurve([R(fill(a)) for a in Emin.a_invariants()])
                    sEred = EllipticCurve([sP.residue_field()(fill(a))
                                           for a in sEmin.a_invariants()])
                    
                    # Checking whether the reduction of the isogeny is
                    # separable
                    phi0num = phi[0].numerator().change_ring(fill)
                    phi0den = phi[0].denominator().change_ring(fill)
                    vphi0num = min(pAdicsK.valuation(cf) for cf in
                                   phi0num.coefficients())
                    vphi0den = min(pAdicsK.valuation(cf) for cf in
                                   phi0den.coefficients())
                    phi0mul = pAdicsK.uniformizer()^(-min(vphi0num, vphi0den))
                    phi0num = phi0num * phi0mul
                    phi0den = phi0den * phi0mul
                    phi0 = (phi0num.change_ring(R) /
                            phi0den.change_ring(R))
                    F = Rx(phi0(x, 0))
                    if F.derivative(x) == 0:
                        result = None # This will cause the function above
                                      # to try other primes or fail
                    else:
                        # Defining variables needed for both p = 2 and p != 2
                        phi1num = phi[1].numerator().change_ring(fill)
                        phi1den = phi[1].denominator().change_ring(fill)
                        vphi1num = min(pAdicsK.valuation(cf) for cf in
                                       phi1num.coefficients())
                        vphi1den = min(pAdicsK.valuation(cf) for cf in
                                       phi1den.coefficients())
                        phi1mul = pAdicsK.uniformizer()^(-min(vphi1num, vphi1den))
                        phi1num = phi1num * phi1mul
                        phi1den = phi1den * phi1mul
                        phi1 = (phi1num.change_ring(R) /
                                phi1den.change_ring(R))
                        H = phi1(x, 0)
                        G = phi1(x, 1) - H
                        p = P.smallest_integer()
                        c1 = (F - x^p).numerator()

                        # Computing number of points in the set
                        # {P : phi P = Frob P}
                        if p == 2:
                            g = Ered.a1()*x + Ered.a3()
                            h = (x^3 + Ered.a2()*x^2 + Ered.a4()*x +
                                 Ered.a6())
                            c3 = (g*G*h + g*G*H + G^2*h + g^2*H + h^2 +
                                  H^2).numerator()
                            c4 = (g - G).numerator()
                            gc13 = gcd(c1, c3)
                            gc134 = gcd(gc13, c4)
                            gc134g = gcd(gc134, g)
                            num = (1 + gc13.radical().degree() +
                                   gc134.radical().degree() -
                                   gc134g.radical().degree())
                        else:
                            R3 = (4*x^3 + Ered.b2()*x^2 + 2*Ered.b4()*x +
                                  Ered.b6())
                            l = _scalar_of_rational_maps(phi0, phi1, Ered,
                                                         sEred)
                            c2 = (l * R3^((p + 1)/2) -
                                  F.derivative(x) * R3).numerator()
                            num = (1 + 2 * gcd(c1, c2).radical().degree() -
                                   gcd(c1, R3).radical().degree())

                        # Computing a_p(E) and the final result
                        apE = F.numerator().degree() + p - num
                        result = beta(Frob)^(-1) * apE

                    # Adding a result to the dictionary with its
                    # corresponding tree
                    if not (result in results):
                        results[result] = pAdicNode(pAdics=node.pAdics(),
                                                    width=node.width)
                    results[result].merge(node)

                # Turning the results dictionary into a value to
                # return
                Tls = [(r, pAdicTree(variables=self.parameters(), root=T))
                       for r, T in results.items()]
                return ConditionalValue([(r, TreeCondition(T)) for r, T in Tls])

            return apply_to_conditional_value(compute_trace2,
                                              sE.minimal_model(P,
                                                               isomorphism=True,
                                                               condition=con1,
                                                               verbose=verbose,
                                                               precision_cap=precision_cap),
                                              use_condition=True,
                                              default_condition=con1)

        return apply_to_conditional_value(compute_trace1,
                                          self.minimal_model(P,
                                                             isomorphism=True,
                                                             condition=condition,
                                                             verbose=verbose,
                                                             precision_cap=precision_cap),
                                          use_condition=True,
                                          default_condition=condition)

    def _trace_of_frob_mult(self, P, beta, condition, verbose,
                            precision_cap):
        r"""Implementation of :meth:`trace_of_frobenius` in case of
        multiplicative_reduction"""
        K = self.definition_field()
        G = K.galois_group()
        Frob = G.artin_symbol(P)
        if self.degree_map(Frob) != 1:
            return None # The degree of the isogeny should be square
                        # free, hence 1
        p = P.smallest_integer()
        long_text = "The condition that a" + str(p) + "E is "
        short_text = "a" + str(p) + "E == "
        result = beta(Frob)^(-1) * QQ(1 + p)
        return ConditionalValue([(result, condition &
                                  TextCondition(long_text + "+1",
                                                short_text + "+1")),
                                 (-result, condition &
                                  TextCondition(long_text + "-1",
                                                short_text + "-1"))])
            
    def _repr_(self):
        """Give a string representation of a Frey Q-curve.

        .. NOTE::

        This is a direct copy from the code included
        in :class:`EllipticCurve_number_field`

        """
        b = self.ainvs()
        a = [z._coeff_repr() for z in b]
        s = "Frey Q-curve defined by "
        s += "y^2 "
        if a[0] == "-1":
            s += "- x*y "
        elif a[0] == '1':
            s += "+ x*y "
        elif b[0]:
            s += "+ %s*x*y "%a[0]
        if a[2] == "-1":
            s += "- y "
        elif a[2] == '1':
            s += "+ y "
        elif b[2]:
            s += "+ %s*y "%a[2]
        s += "= x^3 "
        if a[1] == "-1":
            s += "- x^2 "
        elif a[1] == '1':
            s += "+ x^2 "
        elif b[1]:
            s += "+ %s*x^2 "%a[1]
        if a[3] == "-1":
            s += "- x "
        elif a[3] == '1':
            s += "+ x "
        elif b[3]:
            s += "+ %s*x "%a[3]
        if a[4] == '-1':
            s += "- 1 "
        elif a[4] == '1':
            s += "+ 1 "
        elif b[4]:
            s += "+ %s "%a[4]
        s = s.replace("+ -","- ")
        s += "over %s "%(self.definition_field(),)
        s += "with parameters %s"%(self.parameters(),)
        return s
