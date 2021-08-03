r"""A class for working with $\Q$-curves

A $\Q$-curve is an elliptic curve defined over some number field that
is isogenous to all its galois conjugates. This file contains the
class Qcurve that represents $\Q$-curves that do not have complex
multiplication.

EXAMPLES::

    sage: from modular_method.elliptic_curves.Qcurves import Qcurve
    sage: K.<t> = QuadraticField(-2)
    sage: R.<x> = K[]
    sage: G.<s> = K.galois_group()
    sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : ((-x^2 - 12*x - 18*(t + 1))/(2*x), t)})
    Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2 with t = 1.414213562373095?*I

Q-curves without CM are modular and are linked to classical newforms
that can be computed::

    sage: from modular_method.elliptic_curves.Qcurves import Qcurve
    sage: K.<t> = QuadraticField(-3)
    sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
    sage: E2 = E.decomposable_twist()
    sage: E2.newform() # long time (12 seconds)
    (q + a0*q^3 + (-a0^2 + 1)*q^5 + O(q^6),
     [Dirichlet character modulo 8 of conductor 8 mapping 7 |--> 1, 5 |--> -1,
      Dirichlet character modulo 1 of conductor 1])

AUTHORS:

- Joey van Langen (2019-03-01): initial version

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

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.ell_number_field import EllipticCurve_number_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

from sage.rings.fraction_field import is_FractionField
from sage.rings.number_field.number_field import CyclotomicField, QuadraticField

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.misc.cachefunc import cached_method, cached_function

from sage.all import Integer, ZZ, QQ, copy
from sage.arith.functions import lcm
from sage.functions.other import sqrt
from sage.functions.generalized import sign

from sage.misc.misc_c import prod as product
from sage.misc.mrange import mrange

from sage.arith.misc import hilbert_symbol, euler_phi, gcd, next_prime

from sage.modular.dirichlet import DirichletGroup

from modular_method.number_fields.galois_group import galois_field_change, galois_field_extend, galois_field_restrict
from modular_method.number_fields.field_constructors import composite_field, _concat_maps, _write_as_im_gen_map, fixed_field, write_as_extension, field_with_root
from modular_method.number_fields.dirichlet_characters import dirichlet_fixed_field, dirichlet_to_galois

from modular_method.group_cohomology.calculations import function_with_coboundary, hilbert90

from modular_method.elliptic_curves.twist import twist_elliptic_curve

from modular_method.modular_forms.newform_wrapper import get_newforms

from modular_method.polynomial.symmetric_polynomials import polynomial_to_symmetric

def _rational_maps_of_urst(u, r, s, t):
    r"""Give the rational maps corresponding to a change of Weierstrass
    model"""
    R = u.parent()
    Rxy = PolynomialRing(R, names=["x", "y"]).fraction_field()
    x, y = Rxy.gens()
    F = x - r
    G = y - s*F - t
    return (F/u**Integer(2) , G/u**Integer(3) )

def _rational_maps_of_isomorphism(phi):
    r"""Give the rational maps corresponding to an isomorphism of Elliptic
    curves

    INPUT:
    
    - ``phi`` -- An isomorphism of elliptic curve given by a
      Weierstrass transformation.

    OUTPUT:

    A tuple consisting of two rational functions $F$ and $G$ in the
    variables $x$ and $y$ giving respectively the $x$-coordinate and
    $y$-coordinate maps of the given isomorphism.

    EXAMPLE::

        sage: from modular_method.elliptic_curves.Qcurves import _rational_maps_of_isomorphism
        sage: E = EllipticCurve([2, 3, 4, 5, 6])
        sage: phi = E.isomorphism_to(E.minimal_model())
        sage: _rational_maps_of_isomorphism(phi)
        (x + 1, x + y + 2)

    """
    urst = phi.tuple()
    return _rational_maps_of_urst(*urst)

def _scalar_of_rational_maps(x_map, y_map, dom, codom):
    r"""Return the scalar associated to an isogeny given by rational maps

    For an isogeny $\phi$ between elliptic curves returns the scalar
    $\lambda$ such that ..MATH

       \phi^* \omega = \lambda \omega,

    where $\omega$ is the invariant differential in the corresponding
    elliptic curve.

    If both elliptic curves are defined over a subfield of the complex
    numbers, this scalar is the same as the scalar in the map ..MATH

       z \mapsto \lambda z

    on the complex numbers that induces this isogeny on the
    corresponding quotients.

    INPUT:

    - ``x_map`` -- The x-coordinate map of the isogeny
    
    - ``y_map`` -- The y-coordinate map of the isogeny
    
    - ``dom`` -- The elliptic curve that is the domain of the isogeny

    - ``codom`` -- The elliptic curve that is the codomain of the
      isogeny

    OUTPUT:

    The unique number $\lambda$ such that..MATH

       \phi^* \omega = \lambda \omega,

    where $\omega$ is the invariant differential in the corresponding
    elliptic curve and $\phi$ is the isogeny defined by the given
    maps.

    """
    R = y_map.parent().base_ring()
    x, y = y_map.parent().gens()
    f1 = dom.defining_polynomial()(x, y, Integer(1) )
    f2 = codom.defining_polynomial()(x, y, Integer(1) )
    lfrac = ((x_map.derivative(x) * f1.derivative(y)) /
             f2.derivative(y)(x_map, y_map))
    if lfrac.numerator().monomials() != lfrac.denominator().monomials():
        raise TypeError(str(lfrac) + " is not an algebraic number " +
                        "as it should be")
    m = lfrac.numerator().monomials()[Integer(0) ]
    l = (lfrac.numerator().monomial_coefficient(m) /
         lfrac.denominator().monomial_coefficient(m))
    if lfrac.numerator() != l * lfrac.denominator():
        raise TypeError(str(lfrac) + " is not an algebraic number " +
                        "as it should be")
    return R(l)

def _scalar_of_isogeny(phi):
    r"""Return the scalar associated to an isogeny.

    For an isogeny $\phi$ between elliptic curves returns the scalar
    $\lambda$ such that ..MATH

       \phi^* \omega = \lambda \omega,

    where $\omega$ is the invariant differential in the corresponding
    elliptic curve.

    If both elliptic curves are defined over a subfield of the complex
    numbers, this scalar is the same as the scalar in the map ..MATH

       z \mapsto \lambda z

    on the complex numbers that induces this isogeny on the
    corresponding quotients.

    INPUT:
    
    - ``phi`` -- An isogeny of elliptic curves

    OUTPUT:

    The unique number $\lambda$ such that..MATH

       \phi^* \omega = \lambda \omega,

    where $\omega$ is the invariant differential in the corresponding
    elliptic curve.

    EXAMPLES::

        sage: from modular_method.elliptic_curves.Qcurves import _scalar_of_isogeny
        sage: E = EllipticCurve([0,0,0,1,0])
        sage: P = E.torsion_points()[0]
        sage: phi = E.isogeny(P)
        sage: _scalar_of_isogeny(phi)
        1

    Note that the scalar of an isogeny and its dual multiply to the
    degree of the isogeny::

        sage: from modular_method.elliptic_curves.Qcurves import _scalar_of_isogeny
        sage: E = EllipticCurve("20a4")
        sage: phi = E.isogenies_prime_degree(3)[0]
        sage: _scalar_of_isogeny(phi)
        3
        sage: _scalar_of_isogeny(phi.dual())
        1
        sage: QQ(_scalar_of_isogeny(phi)) * QQ(_scalar_of_isogeny(phi.dual())) == 3
        True

    """
    Fx, Fy = phi.rational_maps()
    return _scalar_of_rational_maps(Fx, Fy, phi.domain(),
                                    phi.codomain())

class Qcurve_base(EllipticCurve_generic):
    r"""A base class for Q-curves

    A Q-curve is an elliptic curve defined over some number field,
    such that all its galois conjugates are isogeneous to the curve
    itself.

    In this class a Q-curve is represented as an elliptic curve E
    defined over a galois number field K, together with for each
    element s of the galois group of K the x- and y-coordinate maps of
    an isogeny from E to s(E), the conjugate of E by s.

    .. NOTE::

    This class is only intended as a parent class to provide
    functionality. Do not use this class directly.

    This class is intended for Q-curves without complex
    multiplication. Although Q-curves with complex multiplication
    might work, the theory behind many of the methods in this class
    was only intended for Q-curves without complex multiplication.

    EXAMPLE::

        sage: from modular_method.elliptic_curves.Qcurves import Qcurve
        sage: K.<t> = QuadraticField(-2)
        sage: R.<x> = K[]
        sage: G.<s> = K.galois_group()
        sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : ((-x^2 - 12*x - 18*(t + 1))/(2*x), t)})
        Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2 with t = 1.414213562373095?*I

    """
    def _is_cached(self, var):
        return hasattr(self, var) and getattr(self, var) != None

    def __init__(self, curve, isogenies={}, guessed_degrees=[], verbose=False):
        r"""Initialize a Q-curve

        Will build all data associated to a Q-curve in the following
        way.
        
        First of all will initialize the curve itself by using either
        a given curve or by data that can be turned into a curve using
        the function EllipticCurve. The resulting curve should be
        defined over some number field and will be redefined over the
        galois closure of this field.
        
        Next it will compute all the galois conjugates of itself with
        respect to the galois group of its base field.

        Third it will fill in data about the isogenies from the curve
        to a Galois conjugate using the provided data.

        Next it will determine the remaining isogenies using the
        guessed degrees and combining previously found isogenies.
        
        If in this way it can not find isogenies from this curve to
        each Galois conjugate, the initialization will produce an
        error and the resulting object should not be used.

        INPUT:
        
         - ``curve`` -- An elliptic curve over some Galois number
           field or any input that would create such a curve when
           passed to the constructor EllipticCurve. It should be
           defined by a Weierstrass equation of the form .. MATH:

              Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

         - ``isogenies`` -- A dictionary (default: {}) with as keys
           elements s of the Galois group of the base field of the
           Q-curve and as values the corresponding isogeny from this
           curve to the Galois conjugate by s. Such isogenies must be
           defined over the field over which the curve is defined or a
           direct extension thereof. The isogeny can be given as
           either an isogeny as a Sage object; a tuple of a rational
           function in $x$ and an algebraic number, that are
           respectively the $x$-coordinate map of the isogeny and the
           induced scalar multiplication on the differentials; or a
           tuple of three rational functions $F$, $G$, $H$ in $x$ such
           that the isogeny is $(x, y) \mapsto (F(x), G(x) y + H(x))$
           outside points mapping to infinity.

        - ``guessed_degrees`` -- A list (default: []) of strictly
           positive integers indicating possible degrees of isogenies
           between this curve and its Galois conjugates.

        - ``verbose`` -- A boolean value or an integer (default:
           False). When set to True or any value larger then zero will
           print comments to stdout about the computations being done
           whilst busy. If set to False or 0 will not print such
           comments. If set to any negative value will also prevent
           the printing of any warnings. A higher value will cause
           more messages to be printed.

        EXAMPLES:

        One can create a Qcurve using an explicit isogeny::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-2)
            sage: G.<s> = K.galois_group()
            sage: E = EllipticCurve([0, 12, 0, 18*(t + 1), 0])
            sage: Es = EllipticCurve([s(a) for a in E.a_invariants()])
            sage: phi = E.isogeny(E([0,0]), codomain=Es)
            sage: Qcurve(E, isogenies={s : phi})
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2 with t = 1.414213562373095?*I

        However, it is sufficient to simply provide the x-coordinate
        map and the scalar of the isogeny::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-2)
            sage: R.<x> = K[]
            sage: G.<s> = K.galois_group()
            sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : ((-x^2 - 12*x - 18*(t + 1))/(2*x), t)})
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2 with t = 1.414213562373095?*I
        
        It is also possible to make the code 'guess' the isogenies by
        giving a suggestion for their degree::
        
            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?

        Note that giving no data with regards to the isogenies will
        result in an error::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: Qcurve([0, 12, 0, 18*(t + 1), 0])
            Traceback (most recent call last):
            ...
            ValueError: There is not sufficient isogeny information to make [0, 12, 0, 18*t + 18, 0] a Q-curve

        TESTS:

        If the isogeny data does not give a valid isogeny an error is raised::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-2)
            sage: R.<x> = K[]
            sage: G.<s> = K.galois_group()
            sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : ((x^2 - 12*x - 18*(t + 1))/(2*x), t)})
            Traceback (most recent call last):
            ...
            ValueError: The given isogeny data ((1/2*x^2 - 6*x - 9*t - 9)/x, t) for the galois conjugate (1,2) does not give a valid isogeny.

        Filling in isogenies by combining other isogenies works correctly::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<w> = QQ[sqrt(2), sqrt(5)].absolute_field()
            sage: t = sqrt(K(2))
            sage: s = sqrt(K(5))
            sage: a4 = -6 * s^2 * t^2 * (5 + 5*s + 10*t + 5*t^2 + 2*s*t)
            sage: a6 = 8 * (s*t)^3 * (1 + t) * (7 + 15*s + 14*t + 7*t^2 + 6*s*t)
            sage: Qcurve([a4, a6], guessed_degrees=[2, 3])
            Q-curve defined by y^2 = x^3 + (-150*w^3+60*w^2+1950*w-1320)*x + (8200/3*w^3-1960*w^2-107000/3*w+30520) over Number Field in w with defining polynomial x^4 - 14*x^2 + 9

        """
        self._init_curve(curve)
        self._init_isogenies()
        for sigma, phi in isogenies.items():
            self._add_isogeny(sigma, phi)
        flag = self._fill_isogenies()
        for d in guessed_degrees:
            self._add_isogenies_of_degree(d, verbose=verbose)
            flag = self._fill_isogenies()
            if flag:
                break
        if not flag:
            raise ValueError("There is not sufficient isogeny information " +
                             "to make " + str(curve) + " a Q-curve")
        # Check that c^2 is the coboundary of the degree map.
        G = self.definition_field().galois_group()
        d = self.degree_map
        c = self.c
        if not all(d(s) * d(t) / d(s*t) == c(s,t)**Integer(2)  for s in G for t in G):
            raise ValueError("The given isogenies are not valid.")

    def definition_field(self):
        r"""Give the field over which this Q-curve is defined.

        .. NOTE::

        This number field is always Galois.

        OUTPUT:

        The number field over which this Q-curve is defined.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.definition_field()
            Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?

        """
        raise NotImplementedError("Method 'definition_field' not " +
                                  "implemented for 'Qcurve_base'")

    def _init_curve(self, curve):
        r"""Initialize the underlying elliptic curve.

        """
        raise NotImplementedError("Method '_init_curve' not " +
                                  "implemented for 'Qcurve_base'")

    def _galois_cache_key(self, sigma):
        r"""Give a cache key for an element of a galois group"""
        return str(sigma), sigma.parent().number_field()

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
        raise NotImplementedError("Method 'galois_conjugate' not " +
                                  "implemented for 'Qcurve_base'")

    # Isogeny related stuff
    def _init_isogenies(self):
        r"""Initialize the isogeny data.

        """
        # The x & y coordinate maps of isogenies
        self._phi_x = dict()
        self._phi_y = dict()
        # Map from the base field of the elliptic curve:
        K = self.definition_field()
        self._to_Kphi = K.hom(K)
        # Initialize the trivial isogeny that is there:
        R = self.base_ring()
        e = K.galois_group().identity()
        Rxy = PolynomialRing(R, names=["x", "y"]).fraction_field()
        x, y = Rxy.gens()
        self._phi_x[e] = Rxy(x)
        self._phi_y[e] = Rxy(y)

    def _check_isogeny(self, dom, codom, phi_x, phi_y):
        r"""Check if phi_x and phi_y define a valid isogeny from dom to codom"""
        Rxy = phi_x.parent().base()
        x, y = Rxy.gens()
        phi_xn = phi_x.numerator()
        phi_xd = phi_x.denominator()
        phi_yn = phi_y.numerator()
        phi_yd = phi_y.denominator()
        f = codom.defining_polynomial()(phi_xn * phi_yd,
                                        phi_xd * phi_yn,
                                        phi_xd * phi_yd)
        cf = sum(f.coefficient(list(e)) * x**e[Integer(0) ]
                 for e in f.exponents()
                 if e[Integer(1) ] == Integer(2) )
        return f == cf * dom.defining_polynomial()(x, y, Integer(1) )
        
    def _add_isogeny(self, sigma, phi):
        r"""Add an isogeny to the stored isogeny data.

        INPUT:

        - ``sigma`` -- A galois homomorphism of the field over which
          this Q-curve is defined.

        - ``phi`` -- An isogeny from this curve to the galois
          conjugate of this curve by sigma, a tuple of the
          corresponding x-coordinate map and a scalar, or a tuple of
          the rational functions in x defining the isogeny.

        """
        if isinstance(phi, tuple):
            R = phi[Integer(0) ].base_ring()
            if is_FractionField(R):
                R = R.base()
            Rxy = PolynomialRing(R, names=['x','y']).fraction_field()
            x, y = Rxy.gens()
            phi_x = phi[Integer(0) ](x)
            dom = self
            codom = self.galois_conjugate(sigma)
            if len(phi) == Integer(2) :
                phi_y = Rxy(((phi_x.derivative(x)
                              * (Integer(2) *y + dom.a1()*x + dom.a3())
                              / R(phi[Integer(1) ]))
                             - (codom.a1()*phi_x + codom.a3())) / Integer(2) )
            elif len(phi) == Integer(3) : 
                phi_y = Rxy(phi[Integer(1) ](x)*y + phi[Integer(2) ](x))
            # Check if these define a valid isogeny before registering
            if not self._check_isogeny(dom, codom, phi_x, phi_y):
                raise ValueError("The given isogeny data " +
                                 str(phi) +
                                 " for the galois conjugate " +
                                 str(sigma) +
                                 " does not give a valid isogeny.")
            self._phi_x[sigma] = phi_x
            self._phi_y[sigma] = phi_y
        else:
            self._add_isogeny(sigma, (phi.x_rational_map(), _scalar_of_isogeny(phi)))

    def _has_isogeny(self, sigma):
        r"""Tell if an isogeny to a Galois conjugate is already registered.

        """
        return (sigma in self._phi_x and
                self._phi_x[sigma] != None and
                sigma in self._phi_y and
                self._phi_y[sigma] != None)

    def _update_isogeny(self, sigma, change):
        r"""Update the field of one specific isogeny"""
        phi_x, phi_y = self._get_isogeny(sigma, change=change)
        self._phi_x[sigma] = phi_x
        self._phi_y[sigma] = phi_y

    def _update_isogeny_field(self):
        r"""Update the field over which all isogenies are defined

        """
        G = list(self.definition_field().galois_group())
        for i in range(len(G)):
            if self._has_isogeny(G[i]):
                Ki = self._get_isogeny_field(G[i])
                Kphi = self.complete_definition_field()
                if Ki != Kphi:
                    # Put old value in second place, so it will not
                    # change if the composite field is isomorphic.
                    data = composite_field(Ki, self._to_Kphi,
                                           give_maps=True,
                                           names=Ki.variable_name())
                    Kphi, Ki_to_new, old_to_new = data
                    if not Kphi.is_absolute():
                        Kphi = Kphi.absolute_field(names=Kphi.variable_name())
                        old_to_new = _concat_maps(old_to_new, Kphi.structure()[Integer(1) ])
                        Ki_to_new = _concat_maps(Ki_to_new, Kphi.structure()[Integer(1) ])
                    self._to_Kphi = _concat_maps(self._to_Kphi, old_to_new)
                    self._update_isogeny(G[i], Ki_to_new)
                    for j in range(i):
                        if self._has_isogeny(G[j]):
                            self._update_isogeny(G[j], old_to_new)
        if not self.complete_definition_field().is_galois():
            Kphi = self.complete_definition_field().galois_closure(names=('al',)); (al,) = Kphi._first_ngens(1)
            clos = self.complete_definition_field().embeddings(Kphi)[Integer(0) ]
            self._to_Kphi = _concat_maps(self._to_Kphi, clos)
            for s in G:
                if self._has_isogeny(G[i]):
                    self._update_isogeny(G[i], clos)

    def _fill_isogenies(self):
        r"""Attempt to fill in missing isogenies by combining known ones.

        """
        self._update_isogeny_field()
        G = self.definition_field().galois_group()
        Kphi = self.complete_definition_field()
        for s in G:
            for t in G:
                if (not self._has_isogeny(s*t) and
                    self._has_isogeny(s) and
                    self._has_isogeny(t)):
                    sL = galois_field_extend(s, Kphi,
                                             embedding=self._to_Kphi)
                    sL_phi_t = self._get_isogeny(t, change=sL.as_hom())
                    phi_s = self._get_isogeny(s)
                    phi_x = sL_phi_t[Integer(0) ](phi_s)
                    phi_y = sL_phi_t[Integer(1) ](phi_s)
                    dom = self.galois_conjugate(G[Integer(0) ], change_ring=self._to_Kphi)
                    codom = self.galois_conjugate(s*t, change_ring=self._to_Kphi)
                    if not self._check_isogeny(dom, codom, phi_x, phi_y):
                        raise ValueError("Concatenated isogeny is not a valid isogeny.")
                    self._phi_x[s*t] = phi_x
                    self._phi_y[s*t] = phi_y
        return all(self._has_isogeny(s) for s in G)

    def _get_isogeny_field(self, sigma):
        r"""Get the field over which an isogeny is defined"""
        raise NotImplementedError("Method '_get_isogeny_field' not " +
                                  "implemented for 'Qcurve_base'")

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
        raise NotImplementedError("Method '_get_isogeny' not " +
                                  "implemented for 'Qcurve_base'")

    def _add_isogenies_of_degree(self, degree, verbose=False):
        r"""Attempt to find isogenies of a given degree"""
        raise NotImplementedError("Method '_add_isogenies_of_degree' not " +
                                  "implemented for 'Qcurve_base'")

    @cached_method(key=_galois_cache_key)
    def isogeny_scalar(self, sigma):
        r"""Return the scalar of the isogeny from this curve to its sigma
        conjugate.

        The scalar of an isogeny is the constant $\lambda$ such that
        the pullback of the invariant differential is $\lambda$ times
        the invariant differential.

        INPUT:
        
        - ``sigma`` -- A galois homomorphism of a number field

        OUTPUT:

        An algebraic integer $\lambda$ such that the isogeny from this
        curve to its `sigma` Galois conjugate pulls back the invariant
        differential to $\lambda$ times the invariant differential of
        this curve.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: G = K.galois_group()
            sage: E.isogeny_scalar(G[0])
            1
            sage: E.isogeny_scalar(G[1])
            1/10*lu^3 + 3/10*lu
            sage: E.isogeny_scalar(G[1]).parent()
            Number Field in lu with defining polynomial x^4 - 2*x^2 + 25

        """
        sigma = galois_field_change(sigma, self.definition_field())
        dom = self.galois_conjugate(sigma**Integer(0) , change_ring=self._to_Kphi)
        codom = self.galois_conjugate(sigma, change_ring=self._to_Kphi)
        return _scalar_of_rational_maps(self._phi_x[sigma],
                                        self._phi_y[sigma], dom,
                                        codom)

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
        raise NotImplementedError("Method 'isogeny_x_map' not " +
                                  "implemented for 'Qcurve_base'")

    def complete_definition_field(self):
        r"""Give the field over which the Q-curve is completely defined.

        OUTPUT:

        A number field over which both this elliptic curve and all
        isogenies from its galois conjugates to itself are defined.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-2)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.complete_definition_field()
            Number Field in t with defining polynomial x^2 + 2 with t = 1.414213562373095?*I

        In general this field is bigger than the definition field::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.definition_field()
            Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?
            sage: E.complete_definition_field()
            Number Field in lu with defining polynomial x^4 - 2*x^2 + 25

        """
        return self._to_Kphi.codomain()

    @cached_method(key=_galois_cache_key)
    def degree_map(self, sigma):
        r"""Give the degree of an isogeny from this curve to a Galois conjugate.

        INPUT:

        - ``sigma`` -- A Galois homomorphism of a number field

        OUTPUT:
        
        The degree of the registered isogeny from this curve to its
        `sigma` Galois conjugate.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(2)
            sage: G.<s> = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.degree_map(s)
            2
        
        The isomorphism from the curve to itself will always have a
        degree that is a square::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: G = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.degree_map(G(1))
            1

        One can also use galois homomorphisms not necessarily defined
        over the definition field::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-1)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: L.<u> = CyclotomicField(12)
            sage: G = L.galois_group()
            sage: [E.degree_map(s) for s in G]
            [1, 2, 1, 2]

        """
        sigma = galois_field_change(sigma, self.definition_field())
        Fx = self._phi_x[sigma].numerator()
        x, y = Fx.parent().gens()
        return ZZ(Fx.degree(x))

    @cached_method
    def degree_map_image(self):
        r"""Give the image of the degree map in $\Q^*/(\Q^*)^2$

        OUTPUT:
        
        A list of squarefree integers such that each value of the
        degree map is a square times such an integer and all integers
        in this list differ a square from a value of the degree map.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: G.<s> = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: [E.degree_map(s) for s in G]
            [1, 2]
            sage: E.degree_map_image()
            [1, 2]

        """
        result = []
        d = self.degree_map
        G = self.definition_field().galois_group()
        for s in G:
            val = d(s).squarefree_part()
            if val not in result:
                result.append(val)
        return result

    def degree_field(self):
        r"""Give the fixed field of the degree map.

        OUTPUT:

        The biggest number field such that for each galois
        homomorphism that acts trivially on this field the degree map
        takes a value in $\Q^2$.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.degree_field()
            Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?

        The degree field is always a subfield of the definition field,
        but can be strictly smaller::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = CyclotomicField(3)
            sage: G.<s> = K.galois_group()
            sage: R.<x> = K[]
            sage: F = (1/4*t*x^4 + 3/2*x^2 - 9/4*t - 9/4)/(x^3 - 2*t*x^2 + (3*t + 3)*x)
            sage: l = -2*t
            sage: E = Qcurve([0, -2*t, 0, 3*t+3, 0], isogenies={s : (F, l)})
            sage: E.degree_field()
            Rational Field

        """
        Kerd = []
        d = self.degree_map
        G = self.definition_field().galois_group()
        for s in G:
            if d(s).is_square():
                Kerd.append(s)
        return fixed_field(Kerd)

    def dual_basis(self, a1=None):
        r"""Give a dual basis for the degree map.

        INPUT:

         - ``a1`` -- Optional parameter (default: None). If set to a
           non-square integer of which the square root is part of the
           degree field, will ensure that this is the first entry of
           the first list returned.

        OUTPUT:
        
        A tuple containing
        
        - A list of squarefree integers such that their square roots
          generate the degree field. This list is of minimal length
          with respect to such lists.
        
        - A list of non-negative integers of the same length as the
          first.

        The degree map of this Q-curve on any galois homomorphism that
        only changes the sign of the square root of i-th entry of the
        first list, is equal to a square times the i-th entry of the
        second list.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<w> = QQ[sqrt(2), sqrt(5)].absolute_field()
            sage: t = sqrt(K(2))
            sage: s = sqrt(K(5))
            sage: a4 = -6 * s^2 * t^2 * (5 + 5*s + 10*t + 5*t^2 + 2*s*t)
            sage: a6 = 8 * (s*t)^3 * (1 + t) * (7 + 15*s + 14*t + 7*t^2 + 6*s*t)
            sage: E = Qcurve([a4, a6], guessed_degrees=[2, 3, 6])
            sage: E.dual_basis()
            ([10, 2], [2, 6])
            sage: [(sigma(t*s)/(t*s), sigma(t)/t, E.degree_map(sigma)) for sigma in K.galois_group()]
            [(1, 1, 1), (-1, 1, 2), (-1, -1, 3), (1, -1, 6)]

        One can also fix the first entry::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<w> = QQ[sqrt(-2), sqrt(-3)].absolute_field()
            sage: t = sqrt(K(-2))
            sage: s = sqrt(K(-3))
            sage: a4 = -6 * s^2 * t^2 * (5 + 5*s + 10*t + 5*t^2 + 2*s*t)
            sage: a6 = 8 * (s*t)^3 * (1 + t) * (7 + 15*s + 14*t + 7*t^2 + 6*s*t)
            sage: E = Qcurve([a4, a6], guessed_degrees=[2, 3, 6])
            sage: E.dual_basis(a1=-2)
            ([-2, -3], [3, 2])
            sage: [(sigma(t)/t, sigma(s)/s, E.degree_map(sigma)) for sigma in K.galois_group()]
            [(1, 1, 1), (-1, 1, 3), (1, -1, 2), (-1, -1, 6)]

        """
        if a1 != None:
            a1 = a1.squarefree_part()
        d = self.degree_map
        Kd = self.degree_field()
        ai = []
        products = [Integer(1) ]
        if a1 != None and Kd(a1).is_square():
            ai.append(a1)
            products.append(a1)
        for tmp in Kd.subfields(degree=Integer(2) ):
            a = tmp[Integer(0) ].discriminant().squarefree_part()
            if a not in products:
                ai.append(a)
                products.extend([(a*b).squarefree_part() for b in products])
        di = [Integer(0) ]*len(ai)
        for sigma in Kd.galois_group():
            ls = [sigma(sqrt(Kd(a)))/sqrt(Kd(a)) for a in ai]
            if sum(ls) == len(ls)-Integer(2) : # Precisely one entry == -1
                for i in range(len(ls)):
                    if ls[i] == -Integer(1) :
                        di[i] = d(sigma)
                        break
        return ai, di

    def _galois_cache_key2(self, sigma, tau):
        r"""Give a cache key for a pair of elements of a galois group"""
        return (self._galois_cache_key(sigma),
                self._galois_cache_key(tau))
    
    @cached_method(key=_galois_cache_key2)
    def c(self, sigma, tau):
        r"""Return the value of the 2-cocycle $c$ associated to this Q-curve.

        For two galois homomorphisms $\sigma$ and $\tau$ of the
        absolute galois group of $\Q$ we have two isogenies from this
        curve conjugated by $\sigma \tau$ to itself, which are the
        isogeny $\phi_{\sigma \tau}$ defining the Q-curve structure
        and the isogeny $\phi_{\sigma} \circ \sigma(\phi(\tau))$ where
        $\phi_{\sigma}$ and $\phi_{\tau}$ are the isogenies from the
        $\sigma$ conjugate and the $\tau$ conjugate of this curve to
        itself respectively. The difference between these two curves
        is a non-zero element c(\sigma, \tau) of $\Q \otimes End(E)$.

        If this curve is non-CM we can identify c(\sigma, \tau) with
        an element of $\Q^*$, hence it defines a map $c: G_\Q^2 \to
        \Q^*$ satisfying .. MATH::

            c(\sigma, \tau) \phi_{\sigma, \tau}
            = \phi_{\sigma} \sigma(\phi_{\tau})

        for all $\sigma, \tau \in G_\Q$, where $G_\Q$ is the absolute
        galois group of $\Q$. Furthermore $c$ is a 2-cocycle.

        In practise the function $c$ can be computed by the fact that
        .. MATH::

            c(\sigma, \tau) = \lambda_\sigma \sigma(\lambda_\tau)
            \lambda_{\sigma \tau}^{-1}

        where $\lambda_\sigma$ is the scalar of the isogeny
        $\phi_\sigma$.

        INPUT:

        - ``sigma`` -- A galois homomorphism over a number field

        - ``tau`` -- A galois homomorphism over a number field

        OUTPUT:

        The value $c(\sigma, \tau)$ as an element of $\Q$, where
        $\sigma$ and $\tau$ are extensions of sigma and tau to
        $\bar{\Q}$ respectively.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: G = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: matrix([[E.c(s, t) for t in G] for s in G])
            [ 1  1]
            [ 1 -2]

        Note that the value of this function always squares to the
        coboundary of the degree map::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-5)
            sage: G = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: d = E.degree_map
            sage: matrix([[E.c(s, t)^2 for t in G] for s in G])
            [1 1]
            [1 4]
            sage: matrix([[d(s)*d(t)/d(s*t) for t in G] for s in G])
            [1 1]
            [1 4]

        """
        l = self.isogeny_scalar
        sigma = galois_field_change(sigma, self.complete_definition_field())
        tau = galois_field_change(tau, self.complete_definition_field())
        return QQ(l(sigma) * sigma(l(tau)) * l(sigma*tau)**(-Integer(1) ))

    def c_pm(self, sigma, tau):
        r"""Return the sign of the 2-cocycle $c$.

        For two galois homomorphisms $\sigma$ and $\tau$ of the
        absolute galois group of $\Q$ we have two isogenies from this
        curve conjugated by $\sigma \tau$ to itself, which are the
        isogeny $\phi_{\sigma \tau}$ defining the Q-curve structure
        and the isogeny $\phi_{\sigma} \circ \sigma(\phi(\tau))$ where
        $\phi_{\sigma}$ and $\phi_{\tau}$ are the isogenies from the
        $\sigma$ conjugate and the $\tau$ conjugate of this curve to
        itself respectively. The difference between these two curves
        is a non-zero element c(\sigma, \tau) of $\Q \otimes End(E)$.

        If this curve is non-CM we can identify c(\sigma, \tau) with
        an element of $\Q^*$, hence it defines a map $c: G_\Q^2 \to
        \Q^*$ satisfying .. MATH::

            c(\sigma, \tau) \phi_{\sigma, \tau}
            = \phi_{\sigma} \sigma(\phi_{\tau})

        for all $\sigma, \tau \in G_\Q$, where $G_\Q$ is the absolute
        galois group of $\Q$. Furthermore $c$ is a 2-cocycle.

        In practise the function $c$ can be computed by the fact that
        .. MATH::

            c(\sigma, \tau) = \lambda_\sigma \sigma(\lambda_\tau)
            \lambda_{\sigma \tau}^{-1}

        where $\lambda_\sigma$ is the scalar of the isogeny
        $\phi_\sigma$.

        INPUT:

        - ``sigma`` -- A galois homomorphism over a number field

        - ``tau`` -- A galois homomorphism over a number field

        OUTPUT:

        The sign of $c(\sigma, \tau)$ as an element of $\Q$, where
        $\sigma$ and $\tau$ are extensions of sigma and tau to
        $\bar{\Q}$ respectively.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: G = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: matrix([[E.c(s, t) for t in G] for s in G])
            [ 1  1]
            [ 1 -2]
            sage: matrix([[E.c_pm(s, t) for t in G] for s in G])
            [ 1  1]
            [ 1 -1]
            sage: matrix([[E.c_abs(s, t) for t in G] for s in G])
            [1 1]
            [1 2]
            sage: all(E.c(s, t) == E.c_pm(s, t) * E.c_abs(s, t) for s in G for t in G)
            True

        """
        return sign(self.c(sigma,tau))

    def c_abs(self, sigma, tau):
        r"""Return the absolute value of the 2-cocycle $c$.

        For two galois homomorphisms $\sigma$ and $\tau$ of the
        absolute galois group of $\Q$ we have two isogenies from this
        curve conjugated by $\sigma \tau$ to itself, which are the
        isogeny $\phi_{\sigma \tau}$ defining the Q-curve structure
        and the isogeny $\phi_{\sigma} \circ \sigma(\phi(\tau))$ where
        $\phi_{\sigma}$ and $\phi_{\tau}$ are the isogenies from the
        $\sigma$ conjugate and the $\tau$ conjugate of this curve to
        itself respectively. The difference between these two curves
        is a non-zero element c(\sigma, \tau) of $\Q \otimes End(E)$.

        If this curve is non-CM we can identify c(\sigma, \tau) with
        an element of $\Q^*$, hence it defines a map $c: G_\Q^2 \to
        \Q^*$ satisfying .. MATH::

            c(\sigma, \tau) \phi_{\sigma, \tau}
            = \phi_{\sigma} \sigma(\phi_{\tau})

        for all $\sigma, \tau \in G_\Q$, where $G_\Q$ is the absolute
        galois group of $\Q$. Furthermore $c$ is a 2-cocycle.

        In practise the function $c$ can be computed by the fact that
        .. MATH::

            c(\sigma, \tau) = \lambda_\sigma \sigma(\lambda_\tau)
            \lambda_{\sigma \tau}^{-1}

        where $\lambda_\sigma$ is the scalar of the isogeny
        $\phi_\sigma$.

        INPUT:

        - ``sigma`` -- A galois homomorphism over a number field

        - ``tau`` -- A galois homomorphism over a number field

        OUTPUT:

        The absolute value of $c(\sigma, \tau)$ as an element of $\Q$,
        where $\sigma$ and $\tau$ are extensions of sigma and tau to
        $\bar{\Q}$ respectively.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: G = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: matrix([[E.c(s, t) for t in G] for s in G])
            [ 1  1]
            [ 1 -2]
            sage: matrix([[E.c_pm(s, t) for t in G] for s in G])
            [ 1  1]
            [ 1 -1]
            sage: matrix([[E.c_abs(s, t) for t in G] for s in G])
            [1 1]
            [1 2]
            sage: all(E.c(s, t) == E.c_pm(s, t) * E.c_abs(s, t) for s in G for t in G)
            True

        """
        return abs(self.c(sigma,tau))

    @cached_method
    def xi_pm(self):
        r"""Return the brauer group representation of the invariant $\xi_\pm$.

        The element $\xi_\pm \in Br_2(\Q)$ is the element in the
        Brauer group corresponding to the homology class of $c_\pm$,
        the cocyle returned by the method :meth:c_pm

        OUTPUT:

        A list of tuples of integers, such that $\xi_\pm$ as an
        element of $Br_2(\Q)$ is the product of the quaternion
        algebras given by each of these tuples.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(2)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.xi_pm()
            [(2, 2)]

        """
        ai, di = self.dual_basis()
        return [(ai[i], di[i]) for i in range(len(ai))]

    @cached_method
    def _xi_pm_primes(self):
        r"""Give the primes at which the $\xi_\pm$ might locally not be 1

        The element $\xi_\pm \in Br_2(\Q)$ is the element in the
        Brauer group corresponding to the homology class of $c_\pm$,
        the cocyle returned by the method :meth:c_pm
        
        OUTPUT:

        A list of prime numbers, such that all primes at which xi_pm
        is locally not trivial are contained in this list. This list
        may contain more primes as well, but not less.

        """
        result = lcm([lcm(h) for h in self.xi_pm()]).prime_factors()
        if Integer(2)  not in result:
            result.insert(Integer(0) ,Integer(2) )
        return result

    @cached_method
    def xi_pm_local(self, p):
        r"""Give a representative for $\xi_\pm$ in $Br_2(\Q_p)$.

        The element $\xi_\pm \in Br_2(\Q)$ is the element in the
        Brauer group corresponding to the homology class of $c_\pm$,
        the cocyle returned by the method :meth:c_pm

        INPUT:

        - ``p`` -- A prime number.

        OUTPUT:
        
        +1 or -1 depending on whether the central simple algebra
        associated to $\xi_\pm$ over $\Q_p$ is split or non-split
        respectively.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.xi_pm()
            [(3, 2)]
            sage: E.xi_pm_local(2)
            -1
            sage: E.xi_pm_local(3)
            -1
            sage: E.xi_pm_local(5)
            1

        """
        if p not in self._xi_pm_primes():
            return Integer(1) 
        else:
            return product([hilbert_symbol(ai,di,p)
                            for (ai,di) in self.xi_pm()])

    def _first_splitting_character(self):
        r"""Compute the first splitting character"""
        N = Integer(1) 
        primes = [p for p in self._xi_pm_primes() if self.xi_pm_local(p) == -Integer(1) ]
        L = CyclotomicField(lcm(euler_phi(p) for p in primes))
        eps_ls = [DirichletGroup(Integer(1) , base_ring=L)[Integer(0) ]]
        for p in primes:
            if p == Integer(2) :
                N *= Integer(4) 
                eps_ls.append(DirichletGroup(Integer(4) , base_ring=L).gen())
            else:
                N *= p
                eps_ls.append(DirichletGroup(p, base_ring=L).gen()**((p-Integer(1) ).odd_part()))
        eps = product([eps_p.extend(N) for eps_p in eps_ls])
        eps = eps.primitive_character()
        eps = eps.minimize_base_ring()
        return eps

    def _splitting_character_data(self, i, j):
        r"""Give data related to splitting characters, i.e. for each splitting
        character we have a list containing the following entries

         0 - the splitting character as a dirichlet character

         1 - the fixed field of that splitting character

         2 - the splitting character as a galois character on its
             fixed field

        INPUT:

        - ``i`` -- The index of the splitting character for which data
          should be retrieved. This may also be a list of such indices
          or one of the special keywords 'all', for all splitting
          characters, or 'conjugacy', for all splitting characters up
          to conjugacy.

        - ``j`` -- The index of the data to be retrieved.

        OUTPUT:
        
        The specified data for each index given in i, formatted
        according to how i was formatted.

        """
        if not self._is_cached('_eps') or Integer(0)  not in self._eps:
            self._eps = dict()
            self._eps[Integer(0) ] = [self._first_splitting_character()]
        if not isinstance(i, str) and hasattr(i, "__iter__"):
            return tuple(self._splitting_character_data(ii, j) for ii in i)
        if i in ZZ:
            if i not in self._eps:
                eps0 = self._eps[Integer(0) ][Integer(0) ]
                chi = self.twist_character(i)
                N = lcm(eps0.conductor(),chi.conductor())
                L, i1, i2 = composite_field(eps0.base_ring(),
                                            chi.base_ring(),
                                            give_maps=True)
                D = DirichletGroup(N, base_ring=L)
                epsi = D(eps0.change_ring(i1)) * D(chi.change_ring(i2))**Integer(2) 
                epsi = epsi.primitive_character()
                epsi = epsi.minimize_base_ring()
                self._eps[i] = [epsi]
            if j >= Integer(1)  and len(self._eps[i]) < Integer(2) :
                self._eps[i].append(dirichlet_fixed_field(self._eps[i][Integer(0) ]))
            if j >= Integer(2)  and len(self._eps[i]) < Integer(3) :
                self._eps[i].append(dirichlet_to_galois(self._eps[i][Integer(0) ]))
            return self._eps[i][j]
        if i == 'all':
            return tuple(self._splitting_character_data(ii, j)
                         for ii in range(self.number_of_splitting_maps()))
        if i == 'conjugacy':
            return tuple(self._splitting_character_data(ii[Integer(0) ], j)
                         for ii in self._conjugacy_classes())
        raise Exception("Invalid index %s."%i)
    
    def splitting_character(self, index=Integer(0) , galois=False):
        r"""Give a splitting character of this Q-curve.

        A splitting character is the character $\eps$ associated to a
        splitting map $\beta$ by $\eps \cdot d = \beta^2$, where $d$
        is the degree map of this Q-curve.

        Although the definition using a splitting map defines a
        splitting character as a galois character, since its fixed
        field is abelian, it can be interpreted as a dirichlet
        character by looking at it as a galois character on
        $\Q(\zeta_n)$. This is the default way of presenting a
        splitting character.

        ALGORITHM:

        To compute the first splitting character, we use the relation
        between splitting characters and the invariant $\xi_\pm$
        computed by :meth:`xi_pm`. The other splitting characters are
        computed using the twist characters computed by
        :meth:`twist_character`.

        .. SEEALSO::

            :meth:`degree_map`,
            :meth:`splitting_map`

        INPUT:

        - ``index`` -- The index (default: 0) of the corresponding
          splitting map. Accepted values are non-negative integers
          smaller than the total amount of splitting maps or one of
          the special values: 'all' for a tuple of all splitting
          characters; and 'conjugacy' for a tuple of splitting
          characters, one for each conjugacy class of splitting
          maps. Also accepts tuples of accepted values including
          tuples themselves.

        - ``galois`` -- A boolean (default: False) indicating whether
          the splitting characters should be given as galois or
          dirichlet characters.

        OUTPUT:

        The splitting character of the given index, given as a
        galois character if galois is set to True or as a Dirichlet
        character otherwise. If the index was 'all' or 'conjugacy'
        will return a tuple of such characters corresponding to the
        corresponding tuple of indices. If the given index was a tuple
        will return a tuple of outputs on each entry of this tuple in
        the same order.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.splitting_character()
            Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1

        There are as many splitting characters as the degree of the
        decomposition field::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.splitting_character('all')
            (Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1,
             Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1,
             Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1,
             Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1,
             Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1,
             Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1,
             Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1,
             Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1)
            sage: E.decomposition_field().degree()
            8

        The galois version of the splitting character relates to the
        degree map and the corresponding splitting map::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(12)
            sage: E = Qcurve([0, 6*t - 12, 0, -54*t + 180, 0], guessed_degrees=[2])
            sage: n = E.decomposition_field().degree()
            sage: G = E.decomposition_field().galois_group()
            sage: d = E.degree_map
            sage: eps = E.splitting_character
            sage: beta = E.splitting_map
            sage: all(eps(i, galois=True)(s) * d(s) == beta(i)(s)^2 for s in G for i in range(n))
            True

        """
        if galois:
            return self._splitting_character_data(index, Integer(2) )
        else:
            return self._splitting_character_data(index, Integer(0) )

    def splitting_character_field(self, index=Integer(0) ):
        r"""Give the fixed field of a splitting character of this Q-curve.

        .. SEEALSO::

            :meth:`splitting_character`

        INPUT:

        - ``index`` -- The index (default: 0) of the corresponding
          splitting map. Accepted values are non-negative integers
          smaller than the total amount of splitting maps or one of
          the special values: 'all' for a tuple of all splitting
          character fields; and 'conjugacy' for a tuple of splitting
          character fields, one for each conjugacy class of splitting
          maps. Also accepts tuples of accepted values including
          tuples themselves.

        OUTPUT:

        The fixed field of a splitting character of the given index.
        If the index was 'all' or 'conjugacy' will return a tuple of
        such fields corresponding to the corresponding tuple of
        indices.  If the given index was a tuple will return a tuple
        of outputs on each entry of this tuple in the same order.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.splitting_character_field()
            Number Field in zeta0 with defining polynomial x^2 - 3 with zeta0 = 1.732050807568878?

        In general it is distinct from the complete definition field::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: K1 = E.complete_definition_field(); K1
            Number Field in lu with defining polynomial x^4 - 6*x^2 + 49
            sage: K2 = E.splitting_character_field(); K2
            Number Field in zeta0 with defining polynomial x^4 - 5*x^2 + 5 with zeta0 = 1.902113032590308?
            sage: K1.is_isomorphic(K2)
            False

        """
        return self._splitting_character_data(index, Integer(1) )

    def _splitting_image_field(self, eps, Keps):
        r"""Compute the image field of a splitting map.
        
        INPUT:

        - ``eps`` -- The corresponding dirichlet character.

        - ``Keps`` -- The fixed field of eps as a galois character.

        OUTPUT:
        
        The field in which the corresponding splitting map takes
        values.

        """
        if isinstance(eps, tuple):
            return tuple(self._splitting_image_field(eps[i], Keps[i])
                         for i in range(len(eps)))
        b = None
        if Integer(2) .divides(Keps.degree()):
            b = Keps.subfields(degree=Integer(2) )[Integer(0) ][Integer(0) ].discriminant().squarefree_part()
        ai, di = self.dual_basis(a1=b)
        L = CyclotomicField(Integer(2) *eps.order())
        for i in range(len(di)):
            Kdi = QuadraticField(di[i])
            if ai[i] == b:
                Lbig, L_to_Lbig, Kdi_to_Lbig = composite_field(L, Kdi,
                                                               give_maps=True)
                alpha = L_to_Lbig(L.gen()) * Kdi_to_Lbig(Kdi.gen())
                L = Lbig.subfield(alpha)[Integer(0) ]
            else:
                L = composite_field(L, Kdi)
        return L

    @cached_method
    def splitting_image_field(self, index=Integer(0) ):
        r"""Give the image field of a splitting map of this Q-curve.

        ALGORITHM:

        The field in which a splitting map takes values is computed
        using only the splitting character and its fixed field.

        .. SEEALSO::

            :meth:`splitting_map`,
            :meth:`splitting_character`,
            :meth:`splitting_character_field`

        INPUT:

        - ``index`` -- The index (default: 0) of the corresponding
          splitting map. Accepted values are non-negative integers
          smaller than the total amount of splitting maps or one of
          the special values: 'all' for a tuple of all image fields;
          and 'conjugacy' for a tuple of image fields, one for each
          conjugacy class of splitting maps. Also accepts tuples of
          accepted values including tuples themselves.

        OUTPUT:

        The image field of the splitting map of the given index.  If
        the index was 'all' or 'conjugacy' will return a tuple of such
        fields corresponding to the corresponding tuple of indices. If
        the given index was a tuple will return a tuple of outputs on
        each entry of this tuple in the same order.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12*t - 12, 0, 180 - 108*t, 0], guessed_degrees=[2])
            sage: K1 = E.splitting_image_field(); K1
            Number Field in zeta4a0 with defining polynomial x^2 + 2 with zeta4a0 = -1/2*zeta4a^2 + 1/2
            sage: beta = E.splitting_map()
            sage: G = E.splitting_field().galois_group()
            sage: [beta(s) for s in G]
            [1, zeta4a0]

        """
        eps = self.splitting_character(index)
        Keps = self.splitting_character_field(index)
        return self._splitting_image_field(eps, Keps)

    def _splitting_field(self, Keps, names=None):
        r"""Imlementation of :meth:`splitting_field`"""
        if isinstance(Keps, tuple):
            return tuple(self._splitting_field(Keps_i) for Keps_i in Keps)
        Kd = self.degree_field()
        return composite_field(Kd, Keps, names=names)

    @cached_method(key=lambda self, index, names: index)
    def splitting_field(self, index=Integer(0) , names=None):
        r"""Give a splitting field of this Q-curve.

        A splitting map $\beta : G_\Q \to \overline{\Q}$ where $G_\Q$
        is the absolute galois group of $\Q$ and $\overline{\Q}$ is
        the algebraic closure of $\Q$ can be seen as a map to
        $\overline{\Q} / \Q$. Since its coboundary with respect to the
        trivial action of $G_\Q$ takes values in $\Q$, regarding
        $\beta$ in this way makes it a 1-cocycle hence a
        homomorphism. Let $H$ be the kernel of the homomorphism $\beta
        : G_\Q \to \overline{\Q} / \Q$. We call the fixed field of $H$
        the fixed field of $\beta$.

        ALGORITHM:
        
        A splitting field is computed using the fixed field of the
        degree map and the fixed field of the corresponding splitting
        character.

        .. SEEALSO::

            :meth:`degree_map`,
            :meth:`splitting_map`,
            :meth:`splitting_character`,
            :meth:`splitting_character_field`,
            :meth:`degree_field`

        INPUT:

        - ``index`` -- The index (default: 0) of the corresponding
          splitting map. Accepted values are non-negative integers
          smaller than the total amount of splitting maps or one of
          the special values: 'all' for a tuple of all splitting
          fields; and 'conjugacy' for a tuple of splitting fields, one
          for each conjugacy class of splitting maps. Also accepts
          tuples of accepted values including tuples themselves.

        - ``names`` -- A string providing a name for the generator of
          the returned field. Note that this argument is ignored if a
          splitting field has been computed before. If set to None
          will use the default naming provided by
          :func:`composite_field`.

        OUTPUT:

        The splitting field corresponding to the splitting map of the
        given index. If the index was 'all' or 'conjugacy' will return
        a tuple of such fields corresponding to the corresponding
        tuple of indices. If the given index was a tuple will return a
        tuple of outputs on each entry of this tuple in the same
        order.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.splitting_field()
            Number Field in zeta0 with defining polynomial x^4 - 5*x^2 + 5 with zeta0 = 1.902113032590308?
            sage: E.degree_field()
            Number Field in t with defining polynomial x^2 - 5 with t = 2.236067977499790?
            sage: E.splitting_character_field()
            Number Field in zeta0 with defining polynomial x^4 - 5*x^2 + 5 with zeta0 = 1.902113032590308?

        """
        Keps = self.splitting_character_field(index)
        return self._splitting_field(Keps, names=names)

    @cached_method(key=lambda self, names: ())
    def decomposition_field(self, names=None):
        r"""Give the field over which the restriction of scalars of this
        Q-curve could decompose as a product of abelian varieties of
        GL_2-type.

        .. NOTE::

        In order for the restriction of scalars to decompose into a
        product of abelian varieties of GL_2-type over this field, the
        cocycle associated to the Q-curve and the coboundary obtained
        from each splitting map must agree over this field.  By
        definition of the splitting maps the coboundaries may however
        differ by the coboundary of a function from the absolute
        galois group of $\QQ$ with values in $\QQ^*$, which does not
        necessarily arise from a function on the galois group of the
        decomposition field.

        .. SEEALSO::

            :meth:`splitting_map`,
            :meth:`splitting_field`,
            :meth:`complete_definition_field`,
            :meth:`does_decompose`,
            :meth:`c`,
            :meth:`decomposable_twist`

        INPUT:

        - `names` -- A string providing a name for the generator of
          the returned field. Note that this argument is ignored if a
          decomposition field has been computed before. If set to None
          will use the default naming provided by
          :func:`composite_field`.

        OUTPUT:
        
        The composite field of :meth:`complete_definition_field` and
        :meth:`splitting_field`

        EXAMPLE::
        
            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.complete_definition_field()
            Number Field in lu with defining polynomial x^4 + 10*x^2 + 1
            sage: E.splitting_field()
            Number Field in tzeta0 with defining polynomial x^4 + 36
            sage: E.decomposition_field()
            Number Field in tzeta0lu with defining polynomial x^8 - x^4 + 1

        """
        Ksplit = self.splitting_field()
        K = self.definition_field()
        Kd = self.degree_field()
        to_Kphi = self._to_Kphi * Kd.embeddings(K)[Integer(0) ]
        to_Ksplit = Kd.embeddings(Ksplit)[Integer(0) ]
        # Put to_Kphi in the second place, so the decomposition field
        # is actually equal to the complete definition field if it
        # would be isomorphic.
        Kdec, from_Ksplit, from_Kphi = composite_field(to_Ksplit,
                                                       to_Kphi,
                                                       give_maps=True,
                                                       names=names)
        self._to_Kdec = from_Kphi
        if not Kdec.is_absolute():
            Kdec = Kdec.absolute_field(names=Kdec.variable_name())
            self._to_Kdec = _concat_maps(from_Kphi, Kdec.structure()[Integer(1) ])
        return Kdec

    @cached_method
    def _over_Kdec(self):
        "Return this Q-curve base changed to the decomposition field."
        if self.definition_field() == self.decomposition_field():
            return self
        else:
            iota = self._to_Kdec * self._to_Kphi
            return self.change_ring(iota)

    @cached_method
    def does_decompose(self):
        r"""Determine whether this Q-curve decomposes over its decomposition
        field.

        The restriction of scalars of this Q-curve over its
        decomposition field decomposes into abelian varieties of
        GL_2-type if and only if the coboundaries of the splitting
        maps differ from the cocycle associated to this Q-curve by a
        coboundary of a function on the absolute galois group of $\QQ$
        with values in $\QQ^*$ that comes from a function defined on
        the galois group of the decomposition field.

        .. NOTE::

        In the case that this Q-curve does decompose, the splitting
        maps will be initialized in such a way that their coboundary
        actually agrees with the cocycle associated to this Q-curve.

        Furthermore, any Q-curve has a twist over its decomposition field
        that does decompose over that field. This can be computed using
        the method :meth:`decomposable_twist`.

        .. SEEALSO::

            :meth:`decomposition_field`,
            :meth:`c`
            :meth:`c_splitting_map`,
            :meth:`splitting_map`,
            :meth:`decomposable_twist`

        OUTPUT:
        
        True if the restriction of scalars of this Q-curve over the
        decomposition field is isogeneous over $\Q$ to the product of
        $\Q$-simple, non-$\Q$-isogeneous abelian varieties of
        $GL_2$-type.

        False otherwise.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.does_decompose()
            False
            sage: E.decomposable_twist().does_decompose()
            True

        """
        if not self._is_cached("_beta"):
            self.splitting_map(verbose=-Integer(1) )
        c = self.c
        c_beta = self.c_splitting_map
        G = self.decomposition_field().galois_group()
        return all(c(s,t) == c_beta(s,t) for s in G for t in G)

    def _splitting_map_first_guess(self):
        r"""Give a naive guess of a splitting map.

        """
        Leps = self.splitting_character().base_ring()
        eps = self.splitting_character(galois=True)
        d = self.degree_map
        Lbeta = self.splitting_image_field()
        iota = Leps.embeddings(Lbeta)[Integer(0) ]
        @cached_function(key=lambda s: (str(s), s.parent().number_field()))
        def beta(sigma):
            return sqrt(iota(d(sigma) * eps(sigma)))
        return beta

    def _first_splitting_map(self, verbose=False):
        r"""Compute a splitting map corresponding to the first splitting
        character.

        """
        self._beta = self._splitting_map_first_guess()
        G = self.decomposition_field().galois_group()
        def c_err(sigma, tau):
            return QQ(self.c(sigma, tau) / self.c_splitting_map(sigma, tau))
        def convert(a):
            if a == Integer(1) :
                return [Integer(0) ]
            elif a == -Integer(1) :
                return [Integer(1) ]
            else:
                raise ValueError("%s is not 1 or -1"%a)
        try:
            alpha = function_with_coboundary(G, (Integer(1) , [-Integer(1) ], [Integer(2) ], convert), c_err)
            beta0 = self._beta
            @cached_function(key=lambda s: (str(s), s.parent().number_field()))
            def beta(sigma):
                sigma = galois_field_change(sigma, self.decomposition_field())
                return beta0(sigma) * alpha(sigma)
            self._beta = beta
            self.c_splitting_map.clear_cache() # Delete values of previous beta
            if not self.does_decompose():
                raise ValueError("Should be impossible to reach this code!");
        except ArithmeticError:
            if verbose >= Integer(0) :
                print("Warning: The restriction of scalars of this Q-curve " +
                       "over the decomposition field does not decompose into "+
                       "abelian varieties of GL_2-type. Use the method " +
                       "decomposable_twist to find a twist that does.")
        return self._beta

    def _indexed_splitting_map(self, i):
        r"""Compute the i-th splitting map from the first splitting map and the
        corresponding twist.

        """
        beta0 = self.splitting_map()
        Lbeta0 = self.splitting_image_field()
        chi = self.twist_character(i, galois=True)
        Lchi = self.twist_character(i).base_ring()
        Lbeta = self.splitting_image_field(i)
        L, from_Lbeta0, from_Lchi = composite_field(Lbeta0, Lchi,
                                                    give_maps=True)
        L, to_L, from_Lbeta = composite_field(L, Lbeta, give_maps=True)
        L, conv = write_as_extension(from_Lbeta, give_map=True)
        to_L = conv * to_L
        from_Lbeta0 = to_L * from_Lbeta0
        from_Lchi = to_L * from_Lchi
        @cached_function(key=lambda s: (str(s), s.parent().number_field()))
        def beta(sigma):
            result = (from_Lbeta0(beta0(sigma)) * from_Lchi(chi(sigma))).list()
            if not all(result[i] == Integer(0)  for i in range(Integer(1) , len(result))):
                raise ValueError("Value " + str(L(result)) +
                                 "of splitting map is not an" +
                                 "element of the image field" +
                                 str(Lbeta))
            return result[Integer(0) ]
        return beta

    @cached_method(key=lambda self, i, v: i)
    def splitting_map(self, index=Integer(0) , verbose=False):
        r"""Give a splitting map of this Q-curve.

        A splitting map is a function from the absolute galois group
        of $\QQ$ to the units in the algebraic closure of $\QQ$ such
        that its coboundary with respect to the trivial action of the
        absolute galois group gives the same cohomology class as the
        cocycle associated to this Q-curve in $H^2(G_\QQ,
        \QQ^*)$. Here $G_\QQ$ is the absolute galois group of $\Q$.

        This method only returns splitting maps that arise as
        functions on the galois group of the decomposition
        field. Furthermore all splitting maps returned by this
        function will have the same coboundary.

        .. NOTE::

        In the case that this Q-curve decomposes over its
        decomposition field the splitting maps will be chosen in such
        a way that the common coboundary is precisely the cocycle
        associated to this Q-curve.

        .. SEEALSO::

            :meth:`splitting_character`,
            :meth:`c`
            :meth:`c_splitting_map`,
            :meth:`does_decompose`,
            :meth:`decomposable_twist`

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
          map. Accepted values are non-negative integers smaller than
          the total amount of splitting maps or one of the special
          values: 'all' for a tuple of all splitting maps; and
          'conjugacy' for a tuple of splitting maps, one for each
          conjugacy class of splitting maps. Also accepts tuples of
          accepted values including tuples themselves.

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

        OUTPUT:

        The splitting map of the given index. If the index was 'all'
        or 'conjugacy' will return a tuple of such maps corresponding
        to the corresponding tuple of indices. If the given index was
        a tuple will return a tuple of outputs on each entry of this
        tuple in the same order.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12*t - 12, 0, 180 - 108*t, 0], guessed_degrees=[2])
            sage: beta = E.splitting_map()
            sage: G = E.decomposition_field().galois_group()
            sage: [beta(s) for s in G]
            [1, zeta4a0]
            sage: [beta(s).minpoly() for s in G]
            [x - 1, x^2 + 2]

        Note that the cocycle associated to this Q-curve and the
        coboundary of the splitting map agree if this Q-curve
        decomposes::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.does_decompose()
            False
            sage: G = E.decomposition_field().galois_group()
            sage: matrix([[E.c(s, t) for t in G] for s in G])
            [ 1  1  1  1  1  1  1  1]
            [ 1  2  1  2  2  1  2  1]
            [ 1 -1  1 -1 -1  1 -1  1]
            [ 1  2  1  2  2  1  2  1]
            [ 1 -2  1 -2 -2  1 -2  1]
            [ 1  1  1  1  1  1  1  1]
            [ 1 -2  1 -2 -2  1 -2  1]
            [ 1 -1  1 -1 -1  1 -1  1]
            sage: matrix([[E.c_splitting_map(s, t) for t in G] for s in G])
            [ 1  1  1  1  1  1  1  1]
            [ 1  2 -1  2  2 -1  2  1]
            [ 1 -1 -1  1 -1 -1  1  1]
            [ 1  2  1 -2  2  1 -2  1]
            [ 1  2 -1  2  2 -1  2  1]
            [ 1 -1 -1  1 -1 -1  1  1]
            [ 1  2  1 -2  2  1 -2  1]
            [ 1  1  1  1  1  1  1  1]
            sage: E2 = E.decomposable_twist()
            sage: E2.does_decompose()
            True
            sage: G = E2.decomposition_field().galois_group()
            sage: matrix([[E2.c(s, t) for t in G] for s in G])
            [ 1  1  1  1]
            [ 1  2  2 -1]
            [ 1  2 -2  1]
            [ 1 -1  1 -1]
            sage: matrix([[E2.c_splitting_map(s, t) for t in G] for s in G])
            [ 1  1  1  1]
            [ 1  2  2 -1]
            [ 1  2 -2  1]
            [ 1 -1  1 -1]

        """
        if not isinstance(index, str) and hasattr(index, "__iter__"):
            return tuple(self.splitting_map(i, verbose=verbose) for i in index)
        if index in ZZ:
            if index == Integer(0) :
                return self._first_splitting_map(verbose=verbose)
            else:
                return self._indexed_splitting_map(index)
        if index == 'all':
            n = self.number_of_splitting_maps()
            return self.splitting_map(tuple(range(n)), verbose=verbose)
        if index == 'conjugacy':
            return tuple(self.splitting_map(index=ii[Integer(0) ], verbose=verbose)
                         for ii in self._conjugacy_classes())
        raise Exception("Invalid index %s"%index)

    @cached_method
    def _splitting_field_map(self, index=Integer(0) ):
        r"""Determine how the image of a splitting character embeds into the
        image of the corresponding splitting map.

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
          map. Accepted values are non-negative integers smaller than
          the total amount of splitting maps.

        OUTPUT:

        A field homomorphism from the image field of the splitting
        character with index `index` to the image field of the
        splitting map with index `index`. This homomorphism is such
        that the values of the splitting character and the
        corresponding splitting map correspond in the image field
        according to the relation ..MATH

            \beta^2 = \eps d,

        where $\beta$ is the splitting map, $\eps$ is the splitting
        character and $d$ is the degree map.

        """
        Leps = self.splitting_character(index).base_ring()
        Lbeta = self.splitting_image_field(index)
        eps = self.splitting_character(index, galois=True)
        beta = self.splitting_map(index)
        d = self.degree_map
        K = self.splitting_field()
        G = K.galois_group()
        return [iota for iota in Leps.embeddings(Lbeta)
                if all(d(s) * iota(eps(s)) == beta(s)**Integer(2) 
                       for s in G)][Integer(0) ]

    @cached_method(key=_galois_cache_key2)
    def c_splitting_map(self, sigma, tau):
        r"""Evaluate the coboundary of a splitting map of this Q-curve.

        .. NOTE::

        This is independent of the chosen splitting map.

        .. SEEALSO::

            :meth:`splitting_map`
       
        INPUT:
        
        - ``sigma`` -- A galois homomorphism of a number field.

        - ``tau`` -- A galois homomorphism of a number field.

        - ``verbose`` -- A boolean value or an integer (default:
          False). When set to True or any value larger then zero will
          print comments to stdout about the computations being done
          whilst busy. If set to False or 0 will not print such
          comments.  If set to any negative value will also prevent
          the printing of any warnings. Higher values will cause more
          messages to be printed.

        OUTPUT:
        
        The value .. MATH::

            \beta(\sigma) \cdot \beta(\tau) \cdot \beta(\sigma
            \tau)^{-1}

        for $\beta$ a splitting map of this Q-curve and $\sigma$ and
        $\tau$ galois extensions to $\bar{\Q}$ of sigma and tau
        respectively.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.does_decompose()
            False
            sage: G = E.decomposition_field().galois_group()
            sage: matrix([[E.c(s, t) for t in G] for s in G])
            [ 1  1  1  1  1  1  1  1]
            [ 1  2  1  2  2  1  2  1]
            [ 1 -1  1 -1 -1  1 -1  1]
            [ 1  2  1  2  2  1  2  1]
            [ 1 -2  1 -2 -2  1 -2  1]
            [ 1  1  1  1  1  1  1  1]
            [ 1 -2  1 -2 -2  1 -2  1]
            [ 1 -1  1 -1 -1  1 -1  1]
            sage: matrix([[E.c_splitting_map(s, t) for t in G] for s in G])
            [ 1  1  1  1  1  1  1  1]
            [ 1  2 -1  2  2 -1  2  1]
            [ 1 -1 -1  1 -1 -1  1  1]
            [ 1  2  1 -2  2  1 -2  1]
            [ 1  2 -1  2  2 -1  2  1]
            [ 1 -1 -1  1 -1 -1  1  1]
            [ 1  2  1 -2  2  1 -2  1]
            [ 1  1  1  1  1  1  1  1]
            sage: E2 = E.decomposable_twist()
            sage: E2.does_decompose()
            True
            sage: G = E2.decomposition_field().galois_group()
            sage: matrix([[E2.c(s, t) for t in G] for s in G])
            [ 1  1  1  1]
            [ 1  2  2 -1]
            [ 1  2 -2  1]
            [ 1 -1  1 -1]
            sage: matrix([[E2.c_splitting_map(s, t) for t in G] for s in G])
            [ 1  1  1  1]
            [ 1  2  2 -1]
            [ 1  2 -2  1]
            [ 1 -1  1 -1]

        """
        if not self._is_cached('_beta'):
            self.splitting_map();
        return QQ(self._beta(sigma) * self._beta(tau) *
                  self._beta(sigma*tau)**(-Integer(1) ))

    def _Kphi_roots(self):
        r"""Give a basis for the field of complete definition.

        OUTPUT:
        
        A list of squarefree integers such that the field of complete
        definition is $\Q$ adjoint all roots of these
        integers. Furthermore this list has minimal length in this
        regard.

        """
        Kphi = self.complete_definition_field()
        products = [Integer(1) ]
        result = []
        for tmp in Kphi.subfields(degree=Integer(2) ):
            c = tmp[Integer(0) ].discriminant().squarefree_part()
            if c not in products:
                result.append(c)
                products.extend([(c*b).squarefree_part() for b in products])
        if len(products) < Kphi.degree():
            raise ValueError("This Q-curve is not completely defined over a " +
                             "2-...-2 extension. This method only works " +
                             "when it is.")
        return result

    def cyclotomic_order(self):
        r"""Return the smallest $N$ such that $\Q(\zeta_N)$ contains the
        largest abelian subfield of the decomposition field.

        .. SEEALSO::

            :meth:`decomposition_field`

        OUTPUT:
        
        The smallest non-negative integer $N$ such that the
        decomposition field of this Q-curve as given by
        :meth:`decomposition_field` is a subfield of $\Q(\zeta_N)$.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: N = E.cyclotomic_order(); N
            40
            sage: L = E.decomposition_field(); L
            Number Field in zeta0lu with defining polynomial x^8 - 22*x^6 - 20*x^5 + 199*x^4 + 380*x^3 + 882*x^2 + 740*x + 1721
            sage: len(L.gen().minpoly().change_ring(CyclotomicField(N)).roots()) > 0
            True

        Also works when the decomposition field is not abelian::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: R.<x> = QQ[]
            sage: K.<a> = (x^4 - 2*x^2 - 2).splitting_field()
            sage: E = Qcurve([0, 12, 0, 18*(1 + sqrt(K(3))), 0], guessed_degrees=[2])
            sage: E.decomposition_field().is_abelian()
            False
            sage: E.cyclotomic_order()
            24

        """
        if not self._is_cached('_N') or not self._is_cached('_ker'):
            K = self.decomposition_field()
            if not K.is_abelian():
                G = K.galois_group()
                K = fixed_field([s*t*s**(-Integer(1) )*t**(-Integer(1) ) for s in G for t in G])
            N = K.conductor(check_abelian=False)
            self._N = N
            L = CyclotomicField(N, names=('zeta',)); (zeta,) = L._first_ngens(1)
            a = K.gen().minpoly().change_ring(L).roots()[Integer(0) ][Integer(0) ]
            self._ker = [n for n in range(N) if gcd(n, N) == Integer(1)  and
                         a == sum(ai * zeta**(i*n)
                                  for (i,ai) in enumerate(a.list()))]
        return self._N

    def _init_twist_characters(self):
        r"""Compute all the twist characters to start with.

        """
        N = self.cyclotomic_order()
        ker = self._ker
        D = DirichletGroup(N)
        self._chi = []
        for chi in D:
            for x in ker:
                if chi(x) != Integer(1) :
                    break
            else:
                chi = chi.primitive_character()
                chi = chi.minimize_base_ring()
                self._chi.append([chi])

    def _twist_character_data(self, i, j):
        r"""Give data of the twist characters, storing for each splitting map
        the following data

         0 - The dirichlet character which twists the first splitting
             map into this one.

         1 - The same character as a galois character.

        INPUT:

        - ``i`` -- The index of the twist character for which to
          retrieve the data. This may also be a tuple of accepted
          values or one of the special keywords 'all', for all of
          them, or 'conjugacy', for all of them up to conjugacy.

        - ``j`` -- The index of the data to be retrieved.

        OUTPUT:
        
        The requested data for all indices in i.

        """
        if not self._is_cached('_chi'):
            self._init_twist_characters()
        if not isinstance(i, str) and hasattr(i, "__iter__"):
            return [self._twist_character_data(ii, j) for ii in i]
        if i in ZZ:
            if j == Integer(1)  and len(self._chi[i]) < Integer(2) :
                self._chi[i].append(dirichlet_to_galois(self._chi[i][Integer(0) ]))
            return self._chi[i][j]
        if i == 'all':
            return tuple(self._twist_character_data(ii, j)
                         for ii in range(self.number_of_splitting_maps()))
        if i == 'conjugacy':
            return tuple(self._twist_character_data(ii[Integer(0) ], j)
                         for ii in self._conjugacy_classes())
        raise Exception("Invalid index %s."%i)
    
    def twist_character(self, index=Integer(0) , galois=False):
        r"""Give the twist needed to obtain a certain splitting map.

        Since all splitting maps have the same coboundary, each
        splitting map is the product of the first splitting map with
        some galois character.  As the fixed field of these characters
        is abelian they can also be represented by Dirichlet
        characters, which is the default representation used.

        .. SEEALSO::

            :meth:`splitting_map`

        INPUT:

        - ``index`` -- The index (default: 0) of the corresponding
          splitting map. Accepted values are non-negative integers
          smaller than the total amount of splitting maps or one of
          the special values: 'all' for a tuple of all twist
          characters; and 'conjugacy' for a tuple of twist characters,
          one for each conjugacy class of splitting maps. Also accepts
          tuples of accepted values including tuples themselves.

        - ``galois`` -- A boolean (default: False) indicating whether
          the twist character should be given as a galois or a
          dirichlet character.

        OUTPUT:

        The character such that the first splitting map times this
        character gives the splitting map of the given index. This
        twist character is given as a galois character if galois is
        set to True or as a Dirichlet character otherwise. If the
        index was 'all' or 'conjugacy' will return a tuple of such
        characters corresponding to the corresponding tuple of
        indices. If the given index was a tuple will return a tuple of
        outputs on each entry of this tuple in the same order.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12*t - 12, 0, 180 - 108*t, 0], guessed_degrees=[2])
            sage: beta = E.splitting_map('all')
            sage: xi = E.twist_character('all', galois=True)
            sage: G = E.decomposition_field().galois_group()
            sage: all(beta[i](s) == xi[i](s) * beta[0](s) for s in G for i in range(len(beta)))
            True

        """
        if galois:
            return self._twist_character_data(index, Integer(1) )
        else:
            return self._twist_character_data(index, Integer(0) )

    @cached_method
    def number_of_splitting_maps(self, count_conjugates=True):
        r"""Give the number of splitting maps for this Q-curve.

        For counting the number of splitting maps, only those
        splitting maps are counted that have the same coboundary and
        are defined over the galois group of the decomposition field.

        .. SEEALSO::

            :meth:`splitting_map`,
            :meth:`decomposition_field`,
            :meth:`c_splitting_map`,
            :meth:`splitting_character`,
            :meth:`twist_character`

        INPUT:

        - ``count_conjugates`` -- A boolean (default: True) indicating
          whether conjugate splitting maps should be counted
          separately.

        OUTPUT:
        
        Gives the number of distinct splitting maps of this Q-curve
        defined over the galois group of its decomposition field. If
        the flag count_conjugates is set to False, will return the
        number of galois conjugacy classes of such splitting maps
        instead.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-1)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.number_of_splitting_maps()
            4
            sage: len(E.splitting_character('all'))
            4
            sage: len(E.twist_character('all'))
            4

        """
        if count_conjugates:
            if not self._is_cached('_chi'):
                self._init_twist_characters()
            return len(self._chi)
        else:
            return len(self._conjugacy_classes())

    @cached_method
    def _conjugacy_classes(self):
        r"""Give a tuple of indices that contains for each conjugacy class of
        splitting maps the index of exactly one element thereof.

        """
        beta_ls = self.splitting_map("all")
        beta_del = [beta for beta in beta_ls]
        beta_dict = dict()
        for i in range(len(beta_ls)):
            beta_dict[beta_ls[i]] = i
        Kcomp = self.decomposition_field()
        G = Kcomp.galois_group()
        result = []
        while len(beta_del) > Integer(0) :
            beta0 = beta_del[Integer(0) ]
            L0 = self.splitting_image_field(beta_dict[beta0])
            result.append([])
            for beta in beta_del:
                L = self.splitting_image_field(beta_dict[beta])
                M, L0_to_M, L_to_M = composite_field(L0, L, give_maps=True)
                if M.degree() == L0.degree() and M.degree() == L.degree():
                    for tau in M.galois_group():
                        flag = True
                        for sigma in G:
                            if tau(L0_to_M(beta0(sigma)))!=L_to_M(beta(sigma)):
                                flag = False
                                break
                        if flag:
                            result[-Integer(1) ].append(beta_dict[beta])
                            break
            for j in result[-Integer(1) ]:
                beta_del.remove(beta_ls[j])
        return tuple(result)

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
        raise NotImplementedError("Method '_isogeny_data' not " +
                                  "implemented for 'Qcurve_base'")

    @cached_method
    def minimal_definition_field(self, give_map=False, names=None):
        r"""Give the minimal Galois field over which this curve can be defined.
        
        INPUT:

        - ``give_map`` -- A boolean value (default: False) which
          indicates whether the inclusion map from the minimal
          definition map into the definition field should be returned.

        - ``names`` -- A string or None (default) giving the name of
          the variable of the returned field. If set to None this will
          be the name of the variable in the definition field with "0"
          appended.
        
        OUTPUT:

        A Galois subfield of the definition field, such that this
        curve can be defined over that field. If `give_map` is set to
        True will return a tuple with this subfield as the first entry
        and the inclusion map as the second entry.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: iota = K.embeddings(CyclotomicField(12))[0]
            sage: E2 = E.change_ring(iota)
            sage: E2.minimal_definition_field()
            Number Field in zeta120 with defining polynomial x^2 - 3 with zeta120 = 1.732050807568878?
            sage: E2.definition_field()
            Cyclotomic Field of order 12 and degree 4

        """
        K = self.definition_field()
        G = K.galois_group()
        H = [s for s in G
             if self.a_invariants() == self.galois_conjugate(s).a_invariants()]
        Knew = fixed_field(H)
        if names is None:
            names = Knew.variable_name()
        if not Knew.is_galois():
            Knew = Knew.galois_closure(names=names)
        if give_map:
            return Knew, Knew.embeddings(K)[Integer(0) ]
        else:
            return Knew

    def minimal_complete_definition_field(self, give_maps=False,
                                          names=None, Kmin=None):
        r"""Give the minimal field over which this curve can be completely
        defined.
        
        INPUT:

        - ``give_maps`` -- A boolean value (default: False) which
          indicates whether some field homomorphisms should be
          returned as well.

        - ``names`` -- A string or None (default) giving the name of
          the variable of the returned field. If set to None this will
          be the name of the variable in the complete definition field
          with "0" appended.

        - ``Kmin`` -- A tuple consisting of a Galois number field and
          an embedding of that number field into the definition field
          of this curve, or None (default). The number field should be
          the minimal Galois number field over which this curve is
          defined. If set to None this will be initialized by
          :meth:minimal_definition_field which does not allow naming
          the variable.
        
        OUTPUT:

        A subfield of the complete definition field, such that this
        Q-curve with all its isogenies can be defined over that
        field. If `give_map` is set to True will return a tuple with
        this subfield as the first entry. The second entry will be the
        inclusion map from the minimal definition field to that
        field. The third entry will be the inclusion map from that
        field to the complete definition field of this curve. All
        these maps will respect the inclusions of the minimal
        definition field into the complete definition field.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: iota = K.embeddings(CyclotomicField(12))[0]
            sage: E2 = E.change_ring(iota)
            sage: E2.complete_definition_field()
            Number Field in zeta12lu with defining polynomial x^8 - 18*x^6 + 239*x^4 - 1638*x^2 + 6241
            sage: E2.minimal_complete_definition_field()
            Number Field in zeta12lu00 with defining polynomial x^4 - 38*x^2 + 1225 with zeta12lu00 = 3/3871*zeta12lu0^7 - 29/7742*zeta12lu0^5 + 407/7742*zeta12lu0^3 + 15531/7742*zeta12lu0

        """
        K = self.definition_field()
        if Kmin is None:
            Kmin, Kmin_to_K = self.minimal_definition_field(give_map=True)
        else:
            Kmin, Kmin_to_K = Kmin
        Kphi = self.complete_definition_field()
        if names is None:
            names = Kphi.variable_name() + "0"
        Kmin_to_Kphi = _concat_maps(Kmin_to_K, self._to_Kphi)
        Kbig, to_big = Kphi.galois_closure(names=names, map=True)
        a = to_big(Kmin_to_Kphi(Kmin.gen()))
        G = Kbig.galois_group()
        l = {s : to_big(self.isogeny_scalar(s)) for s in G}
        H = [s for s in G
             if (all(s(ls) == ls for ls in l.values()) and
                 s(a) == a)]
        Kphi_min = fixed_field(H, map=give_maps)
        if give_maps:
            Kphi_min, Kphi_min_to_big = Kphi_min
            b = Kphi_min.gen()
            phi_map = [psi for psi in Kphi_min.embeddings(Kphi)
                       if Kphi_min_to_big(b) == to_big(psi(b))][Integer(0) ]
            min_map = [psi for psi in Kmin.embeddings(Kphi_min)
                       if Kphi_min_to_big(psi(Kmin.gen())) == a][Integer(0) ]
            return Kphi_min, min_map, phi_map
        else:
            return Kphi_min

    def _minimal_a_invariants(self, from_min):
        r"""Give the a_invariants over the minimal definition_field

        INPUT:

        - ``from_min`` -- The map from the minimal definition field to
          the definition field.

        """
        raise NotImplementedError("Method '_minimal_a_invariants' not " +
                                  "implemented for 'Qcurve_base'")

    def _minimal_isogeny_x_maps(self, Kmin, from_min, min_map):
        r"""Give the x-coordinate maps of isogenies for :meth:minimize_fields

        INPUT:

        - ``Kmin`` -- The minimal definition field
        
        - ``from_min`` -- The map from the minimal definition field to
          the definition field.

        - ``min_map`` -- The map from the minimal definition field to
          the minimal complete definition field.

        """
        raise NotImplementedError("Method '_minimal_isogeny_x_maps' not " +
                                  "implemented for 'Qcurve_base'")
    
    def _minimal_isogeny_scalars(self, Kmin, phi_map):
        r"""Give the scalars of isogenies for :meth:minimize_fields

        INPUT:

        - ``Kmin`` -- The minimal definition field

        - ``phi_map`` -- The map from the minimal complete definition
          field to the complete definition field.

        """
        G = Kmin.galois_group()
        l = self.isogeny_scalar
        _, iota = write_as_extension(phi_map, give_map=True)
        return {s : iota(l(s)).list()[Integer(0) ] for s in G}

    def _minimize_fields(self, names=None):
        r"""Give the fields required for :meth:minimize_fields"""
        if names is None:
            names = [None, None]
            
        # Find the necessary fields and maps
        Kmin, from_min = self.minimal_definition_field(names=names[Integer(0) ],
                                                       give_map=True)
        Kphi, min_map, phi_map = self.minimal_complete_definition_field(names=names[Integer(1) ],
                                                                        give_maps=True,
                                                                        Kmin=(Kmin,
                                                                              from_min))

        # Write Kphi as an extension of Kmin
        test = Kphi.gen()
        Kphi, conv = write_as_extension(min_map, give_map=True, names=names[Integer(1) ])
        conv2 = [psi for psi in Kphi.embeddings(test.parent())
                 if psi(conv(test)) == test][Integer(0) ]
        min_map = _concat_maps(min_map, conv)
        phi_map = _concat_maps(conv2, phi_map)

        return Kmin, Kphi, from_min, min_map, phi_map
        
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
        raise NotImplementedError("Method 'minimize_fields' not " +
                                  "implemented for 'Qcurve_base'")

    def _twist(self, gamma):
        r"""Helper function for :meth:twist"""
        K_gamma = gamma.parent()
        K, iota, gamma_map = composite_field(self._to_Kphi,
                                             K_gamma,
                                             give_maps=True)
        if not K.is_absolute():
            K = K.absolute_field(names=K.variable_name())
            iota = _concat_maps(iota, K.structure()[Integer(1) ])
            gamma_map = _concat_maps(gamma_map, K.structure()[Integer(1) ])
        gamma = gamma_map(gamma)
        E_map = _concat_maps(self._to_Kphi, iota)
        E0 = self.galois_conjugate(self.definition_field().galois_group().identity(),
                                   change_ring=E_map)
        E = twist_elliptic_curve(E0, gamma)
        ainvs = E.a_invariants()
        G = K.galois_group()
        l = self.isogeny_scalar
        F = self.isogeny_x_map
        isogenies=dict()
        R = K['y']; (y,) = R._first_ngens(1)
        for s in G:
            ls = iota(l(s))
            if (gamma/s(gamma)).is_square():
                agamma = sqrt(gamma/s(gamma))
                Fmap = E_map
            else:
                Ls = K.extension(y**Integer(2)  - gamma/s(gamma), names=('agamma',)); (agamma,) = Ls._first_ngens(1)
                Fmap = _concat_maps(E_map, K.hom(Ls))
            Fs = F(s, change_ring=Fmap)
            x = Fs.parent().gen()
            Fs = s(gamma) * Fs(x / gamma)
            ls = agamma * ls
            isogenies[s] = (Fs, ls)
        return ainvs, isogenies

    def twist(self, gamma):
        r"""Give the twist of this Q-curve by a given element gamma.

        If this Q-curve was given by .. MATH::

            E : y^2 = x^3 + a_2 x^2 + a_4 x + a_6

        the twisted Q-curve is given by .. MATH::
        
            E : y^2 = x^3 + \gamma a_2 x^2 + \gamma^2 a_4 x
                      + \gamma^3 a_6

        INPUT:

        - ``gamma`` -- An element of a number field.

        OUTPUT:
        
        A Q-curve which is the twist of this Q-curve by gamma. This
        curve will be defined over the smallest possible field over
        which it is completely defined.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2]); E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 5 with t = 2.236067977499790?
            sage: E2 = E.twist(t); E2
            Q-curve defined by y^2 = x^3 + (-6*lu0)*x^2 + (-45*lu0+90)*x over Number Field in lu0 with defining polynomial x^2 - 20 with lu0 = -1/7*lu^3 + 13/7*lu
            sage: E2.complete_definition_field()
            Number Field in agamma00 with defining polynomial x^4 - 16*x^3 - 8*x^2 + 576*x - 1264

        """
        raise NotImplementedError("Method 'twist' not " +
                                  "implemented for 'Qcurve_base'")

    def _decomposable_twist_set(self):
        r"""Give the set of ideals required to compute the decomposable twist
        of this curve.

        .. SEE_ALSO::

            :meth:`decomposable_twist`,
            :meth:`does_decompose`,
            :meth:`decomposition_field`,
            :meth:`twist`

        OUTPUT:

        A set of ideals of the decomposition field that is closed
        under the action of the galois group. Furthermore it contains
        generators of the two-torsion of the class group of the
        decomposition field modulo ideals that are the product of all
        prime ideals above a common prime number.

        """
        K = self.decomposition_field()
        CG = K.class_group(proof=False)
        Pgen = [CG(product(K.primes_above(p)))
                for p in K.discriminant().prime_factors()]
        Pord = [P.order() for P in Pgen]
        H = []
        for k in mrange(Pord):
            CI = product(Pgen[i]**k[i] for i in range(len(k)))
            if CI not in H:
                H.append(CI)
        S0 = []
        skip = copy(H)
        for CI in CG:
            if CI not in skip and CI**Integer(2)  in H:
                S0.append(CI.ideal())
                skip.extend([CI * h for h in H])
        G = K.galois_group()
        S = []
        for I in S0:
            for P in I.prime_factors():
                for s in G:
                    sP = s(P)
                    if sP not in S:
                        S.append(sP)
        return S
    
    def _decomposable_twist(self):
        r"""Give a twist of this Q-curve for which the restriction of scalars
        over the decomposition field decomposes as a product of
        abelian varieties of GL_2-type.

        .. SEE_ALSO::

            :meth:`decomposable_twist`,
            :meth:`does_decompose`,
            :meth:`decomposition_field`,
            :meth:`twist`

        OUTPUT:
        
        An element of the decomposition field such that twisting this
        Q-curve by that element gives a Q-curve of which the
        restriction of scalars over the decomposition field is an
        abelian variety isogenous to a product of $\QQ$-simple,
        non-$\QQ$-isogenous abelian varieties of GL_2-type.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(7)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.does_decompose()
            False
            sage: E.decomposable_twist().does_decompose()
            True

        """
        K = self.decomposition_field()
        if self.does_decompose():
            return K(Integer(1) )
        S = self._decomposable_twist_set()
        G = K.galois_group()
        US = K.S_unit_group(proof=False, S=S)
        def c_err(sigma, tau):
            return US(self.c(sigma, tau) / self.c_splitting_map(sigma, tau))
        alpha = function_with_coboundary(G, US, c_err)
        return hilbert90(K, lambda s: alpha(s)**Integer(2) )
    
    def decomposable_twist(self):
        r"""Give a twist of this Q-curve for which the restriction of scalars
        over the decomposition field decomposes as a product of
        abelian varieties of GL_2-type.

        .. SEE_ALSO::

            :meth:`does_decompose`,
            :meth:`decomposition_field`,
            :meth:`twist`

        OUTPUT:
        
        A Q-curve which is a twist of this curve that has the same
        decomposition field. When taking the restriction of scalars of
        this new curve over the decomposition field the resulting
        abelian variety is isogenous to a product of $\QQ$-simple,
        non-$\QQ$-isogenous abelian varieties of GL_2-type.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(7)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.does_decompose()
            False
            sage: E.decomposable_twist().does_decompose()
            True

        """
        if self.does_decompose():
            return self
        return self.twist(self._decomposable_twist())
        # twist will do the minimization of the field!

    def _complete_definition_twist(self, roots):
        r"""Give a twist of this curve completely defined over a given field.

        .. SEEALSO::

            :meth:`degree_map_image`,
            :meth:`twist`

        INPUT:
        
        - ``roots`` -- A list of rational numbers such that up to sign
          they form generators of the image of the degree map in
          $\QQ^* / (\QQ^*)^2$.
        
        OUTPUT:

        An element of the decomposition field such that twisting this
        Q-curve by that element gives a Q-curve that is completely
        defined over the definition field with all the roots of the
        integers in the given list `roots` adjoined.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2]); E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?
            sage: E.degree_map_image()
            [1, 2]
            sage: K1 = E.complete_definition_field(); K1
            Number Field in lu with defining polynomial x^4 - 2*x^2 + 25
            sage: K1(2).is_square()
            False
            sage: E2 = E.complete_definition_twist([2]); E2
            Q-curve defined by y^2 = x^3 + (-6*lu0)*x^2 + (-27*lu0+54)*x over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu
            sage: K2 = E2.complete_definition_field(); K2
            Number Field in agamma00 with defining polynomial x^4 - 88*x^2 + 400
            sage: K2(2).is_square()
            True

        """
        # Calculate all elements generated by absolute roots mod squares
        # and how to obtain them
        roots_image = {Integer(1)  : []}
        for i in range(len(roots)):
            if abs(roots[i]).squarefree_part() not in roots_image:
                for b in list(roots_image):
                    bi = abs(roots[i]*b).squarefree_part()
                    roots_image[bi] = roots_image[b] + [i]

        # Check if roots is valid:
        d_image = self.degree_map_image()
        flag = (len(roots_image) != len(d_image)) # At least enough elements
        for a in roots: # Check each root plus or minus one associated to d
            if flag:
                break
            flag = abs(a) not in d_image
        if flag:
            raise ValueError("The set " + str(roots) +
                             " does not give a valid set of roots")

        # Let's compute the fields and corresponding embeddings
        Kbase = self.definition_field()
        Kold = self.complete_definition_field()
        Kroots = QQ
        for i, a in enumerate(roots):
            Kroots = field_with_root(Kroots, a, names=('sqrt_a'+str(i)))
        base_to_old = self._to_Kphi
        Knew, base_to_new, roots_to_new = composite_field(Kbase, Kroots,
                                                          give_maps=True)
        Kbig, old_to_big, new_to_big = composite_field(base_to_old,
                                                       base_to_new,
                                                       give_maps=True)
        base_to_big = _concat_maps(base_to_new, new_to_big)
        _, phi = write_as_extension(base_to_big, give_map=True)

        # The map we want as scalars for the new curve
        d = self.degree_map
        @cached_function(key=lambda s: (str(s), s.parent().number_field()))
        def mu(s):
            my_iter = (roots[i] for i in roots_image[d(s).squarefree_part()])
            return sqrt(Knew(product(my_iter)))

        # The correction map
        l = self.isogeny_scalar
        def alpha(s):
            return new_to_big(mu(s))**Integer(2)  / old_to_big(l(s))**Integer(2) 

        # The twist parameter
        gamma = phi(hilbert90(Kbig, alpha))
        gammals = gamma.list()
        if not all(gamma[i] == Integer(0)  for i in range(Integer(1) , len(gammals))):
            raise ValueError(f"Sought twist {gamma} is not defined " +
                             "over the definition field")
        return gammals[Integer(0) ]

    def complete_definition_twist(self, roots):
        r"""Give a twist of this curve completely defined over a given field.

        .. SEEALSO::

            :meth:`degree_map_image`,
            :meth:`twist`

        INPUT:
        
        - ``roots`` -- A list of rational numbers such that up to sign
          they form generators of the image of the degree map in
          $\QQ^* / (\QQ^*)^2$.
        
        OUTPUT:

        A Q-curve that is a twist of this curve and is defined over
        the same definition field as this curve. Furthermore it is
        completely defined over the definition field with all the
        roots of the integers in the given list `roots` adjoined.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2]); E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?
            sage: E.degree_map_image()
            [1, 2]
            sage: K1 = E.complete_definition_field(); K1
            Number Field in lu with defining polynomial x^4 - 2*x^2 + 25
            sage: K1(2).is_square()
            False
            sage: E2 = E.complete_definition_twist([2]); E2
            Q-curve defined by y^2 = x^3 + (-6*lu0)*x^2 + (-27*lu0+54)*x over Number Field in lu0 with defining polynomial x^2 - 12 with lu0 = -1/5*lu^3 + 7/5*lu
            sage: K2 = E2.complete_definition_field(); K2
            Number Field in agamma00 with defining polynomial x^4 - 88*x^2 + 400
            sage: K2(2).is_square()
            True

        """
        return self.twist(self._complete_definition_twist(roots))

    def conductor_restriction_of_scalars(self):
        r"""Give the conductor of the restriction of scalars of this Q-curve.

        OUTPUT:

        The conductor of the restriction of scalars of this curve over
        the decomposition field.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.conductor_restriction_of_scalars()
            5566277615616

        .. SEE_ALSO::

            :meth:`decomposition_field`

        """
        raise NotImplementedError("Method 'conductor_restriction_of_scalars' not " +
                                  "implemented for 'Qcurve_base'")

    def newform_levels(self, N=None):
        r"""Compute the levels of newforms that could be associated to this
        Q-curve.

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
        This function computes all the possible levels for the
        newforms under these constraints.

        INPUT:

        - ``N`` -- A positive integer or None (default: None) which is
          the conductor of the restriction of scalars of this Q-curve
          over the decomposition field. This may also be a factor $M$
          of the conductor $N$ that is coprime to $N / M$, in which
          case this code only computes the part of the levels that is
          coprime to $N / M$. If set to None will be initialized as
          the result of :meth:`conductor_restriction_of_scalars`.

        OUTPUT:

        A list of tuples, each tuple representing one of the options
        for the levels of the newforms associated to this Q-curve. The
        $i$-th entry of such a tuple is the level of a newform
        corresponding to the $i$-th conjugacy class of splitting maps,
        as returned by :meth:`splitting_map`.

        If the given `N` was only part of the conductor of the
        restriction of scalars, then the given levels will only be the
        part of the levels that only contain primes dividing `N`.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: E.newform_levels()
            [(1536,)]

        If the restriction of scalars decomposes as a product of
        abelian varieties, then there are multiple levels::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: E.number_of_splitting_maps(count_conjugates=False)
            2
            sage: E.newform_levels() # Inconsistency in .unit_group()
            [(7200, 7200)]

        The levels for each component might be distinct, in which case
        the list may contain multiple options of how the levels are
        distributed among the components::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: E.number_of_splitting_maps(count_conjugates=False)
            2
            sage: E.newform_levels()
            [(96, 192), (192, 96)]

        """
        if N is None:
            N = self.conductor_restriction_of_scalars()
        eps = [c**(-Integer(1) ) for c in self.splitting_character(index='conjugacy')]
        chi = [c**(-Integer(1) ) for c in self.twist_character(index='conjugacy')]
        M = lcm(character.modulus() for character in (eps + chi))
        L = QQ
        for i in range(len(eps)):
            L, i1, i2 = composite_field(L, eps[i].base_ring(),
                                        give_maps=True)
            eps[i] = eps[i].change_ring(i2)
            for j in range(i):
                eps[j] = eps[j].change_ring(i1)
        for i in range(len(chi)):
            L, i1, i2 = composite_field(L, chi[i].base_ring(),
                                        give_maps=True)
            chi[i] = chi[i].change_ring(i2)
            for j in range(i):
                chi[j] = chi[j].change_ring(i1)
            eps = [c.change_ring(i1) for c in eps]
        D = DirichletGroup(M, base_ring=L)
        eps = [D(character) for character in eps]
        chi = [D(character) for character in chi]
        alpha = tuple(eps[i].conductor() for i in range(len(eps)))
        beta = tuple(tuple((chi[j] * chi[i]**(-Integer(1) )).conductor()
                           for j in range(len(eps))) for i in range(len(eps)))
        gamma = tuple(tuple((chi[j] * chi[i]**(-Integer(1) ) * eps[i]).conductor()
                            for j in range(len(eps))) for i in range(len(eps)))
        d = [Kf.degree()
             for Kf in self.splitting_image_field(index='conjugacy')]
        primes = ZZ(N).prime_factors()
        level_dict = {p : self._newform_levels(prime=p, alpha=alpha, beta=beta,
                                               gamma=gamma, d=d, N=N)
                      for p in primes}
        return [tuple(product(primes[i]**level_dict[primes[i]][x[i]][j]
                              for i in range(len(primes))) # over all primes
                      for j in range(len(d))) # over all newforms
                # over all possible tuples of exponents
                for x in mrange([len(level_dict[p]) for p in primes])]
    
    def _newform_levels(self, N, prime, alpha, beta, gamma, d):
        r"""Implementation of :meth:`newform_levels`.

        Compute the possible orders of a specific prime in the levels
        of the newforms.

        INPUT:

        - ``N`` -- The argument `N`
        
        - ``prime`` -- A prime number

        - ``alpha`` -- A tuple of the conductors of the characters
          associated to the newforms.

        - ``beta`` -- A tuple of tuples of non-negative integers
          containing at index i, j the conductor of the twist that
          turns the i-th newform into the j-th newform.

        - ``gamma`` -- A tuple of tuples of non-negative integers
          containing at index i, j the conductor of the product of the
          twist that turns the i-th newform into the j-th newform and
          the character of the i-th newform.

        - ``d`` -- A tuple of the respective degrees of the fields in
          which the newforms have their coefficients.

        """
        alpha = tuple(alpha[i].ord(prime) for i in range(len(alpha)))
        beta = tuple(tuple(beta[i][j].ord(prime) for j in range(len(beta[i])))
                     for i in range(len(beta)))
        gamma = tuple(tuple(gamma[i][j].ord(prime)
                            for j in range(len(gamma[i])))
                      for i in range(len(gamma)))
        N = N.ord(prime)
        # Small cases
        x_max = max(max(max(beta[i][j] + Integer(1) , beta[i][j] + gamma[i][j])
                        for j in range(len(beta[i])))
                    for i in range(len(beta)))
        x_ls = []
        if sum(d) * x_max >= N: # Only small possibilities
            for x in mrange([x_max + Integer(1) ]*len(alpha)): 
                candidate = (sum(d[i] * x[i] for i in range(len(d))) == N)
                for i in range(len(beta)):
                    if not candidate:
                        break
                    candidate = (alpha[i] <= x[i])
                    for j in range(len(beta[i])):
                        if not candidate:
                            break
                        if beta[i][j] == Integer(0) :
                            candidate = (x[j] == x[i])
                        elif (x[i] > max(beta[i][j] + Integer(1) ,
                                         beta[i][j] + gamma[i][j]) or
                              (x[i] < max(beta[i][j] + Integer(1) ,
                                          beta[i][j] + gamma[i][j]) and
                               gamma[i][j] >= Integer(2) ) or
                              (x[i] == Integer(1)  and alpha[i] == Integer(1)  and
                               beta[i][j] == Integer(1)  and gamma[i][j] == Integer(1) )):
                            candidate = (x[j] == max(x[i], beta[i][j] + Integer(1) ,
                                                     beta[i][j] + gamma[i][j]))
                        else:
                            candidate = (x[j] <= max(x[i], beta[i][j] + Integer(1) ,
                                                     beta[i][j] + gamma[i][j]))
                if candidate:
                    x_ls.append(tuple(x))
        elif sum(d).divides(N): # Big possibility
            x = ZZ(N / sum(d))
            candidate = True
            for i in range(len(d)):
                if not candidate:
                    break
                candidate = (alpha[i] <= x)
            if candidate:
                x_ls.append(tuple(x for i in range(len(d))))
        return x_ls

    @cached_method
    def _trace_power_formula(self, power):
        r"""Give the formula to compute the trace of frobenius to a given
        power.

        The trace of the galois representation of this newform at a
        Frobenius element to the power $n$ can be computed from the
        trace and determinant at the frobenius element with this
        formula

        INPUT:

        - ``power`` -- The power $n$ of the frobenius element.

        """
        R = QQ['x, y']; (x, y,) = R._first_ngens(2)
        return polynomial_to_symmetric(x**power + y**power)

class Qcurve(Qcurve_base, EllipticCurve_number_field):
    r"""A Q-curve over some number field

    A Q-curve is an elliptic curve defined over some number field,
    such that all its galois conjugates are isogeneous to the curve
    itself.

    In this class a Q-curve is represented as an elliptic curve E
    defined over a galois number field K, together with for each
    element s of the galois group of K the x- and y-coordinate maps of
    an isogeny from E to s(E), the conjugate of E by s.

    .. NOTE::

    This class is intended for Q-curves without complex
    multiplication. Although Q-curves with complex multiplication
    might work, the theory behind many of the methods in this class
    was only intended for Q-curves without complex multiplication.

    EXAMPLE::

        sage: from modular_method.elliptic_curves.Qcurves import Qcurve
        sage: K.<t> = QuadraticField(-2)
        sage: R.<x> = K[]
        sage: G.<s> = K.galois_group()
        sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : ((-x^2 - 12*x - 18*(t + 1))/(2*x), t)})
        Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2 with t = 1.414213562373095?*I

    TESTS:

    Check that given degrees work correctly if there are two galois
    conjugates that are twists of one another::

        sage: from modular_method.elliptic_curves.Qcurves import Qcurve
        sage: R.<x> = QQ[]
        sage: K.<t> = NumberField(x^4 - 10*x^2 + 1)
        sage: Qcurve([0, -12*t, 0, 6*t^3 + 18*t^2 + 6*t, 0], guessed_degrees=[2])
        Q-curve defined by y^2 = x^3 + (-12*t)*x^2 + (6*t^3+18*t^2+6*t)*x over Number Field in t with defining polynomial x^4 - 10*x^2 + 1

    """
    def definition_field(self):
        r"""Give the field over which this Q-curve is defined.

        .. NOTE::

        This number field is always Galois.

        OUTPUT:

        The number field over which this Q-curve is defined.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.definition_field()
            Number Field in t with defining polynomial x^2 - 3 with t = 1.732050807568878?

        """
        return self.base_ring()
    
    def _init_curve(self, curve):
        r"""Initialize the underlying elliptic curve.

        """
        if not isinstance(curve, EllipticCurve_number_field):
            curve = EllipticCurve(curve)
        if not isinstance(curve, EllipticCurve_number_field):
            raise ValueError("%s can not be a Q-curve"%curve)
        K = curve.base_ring()
        if not K.is_galois():
            raise ValueError("The curve should be defined over a Galois number field")
        else:
            EllipticCurve_number_field.__init__(self, K, curve.a_invariants())

    @cached_method(key=lambda self, sigma, change :
                   (self._galois_cache_key(sigma), change))
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
        return EllipticCurve(self.a_invariants()).change_ring(change)

    # Isogeny related stuff
    def _get_isogeny_field(self, sigma):
        r"""Get the field over which an isogeny is defined"""
        return self._phi_x[sigma].base_ring()

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
            Sxy = Rxy.change_ring(change.codomain()).fraction_field()
            phi_x = (Sxy(phi_x.numerator().change_ring(change)) /
                     Sxy(phi_x.denominator().change_ring(change)))
            phi_y = (Sxy(phi_y.numerator().change_ring(change)) /
                     Sxy(phi_y.denominator().change_ring(change)))
        return phi_x, phi_y

    def _add_isogenies_of_degree(self, degree, verbose=False):
        r"""Attempt to find isogenies of a given degree.

        """
        G = self.definition_field().galois_group()
        fd = self.torsion_polynomial(degree)
        for g,e in fd.factor():
            Kd = self.definition_field().extension(g, names=('l',)); (l,) = Kd._first_ngens(1)
            Ed = EllipticCurve_number_field.base_extend(self, Kd)
            S = Kd['x']; (x,) = S._first_ngens(1)
            psi = Ed.isogeny(x - l)
            E_t = psi.codomain()
            j_t = E_t.j_invariant()
            for s in G:
                if Kd(self.galois_conjugate(s).j_invariant()) == j_t:
                    if verbose > Integer(0) :
                        print("Degree %s isogeny found for"%degree, s)
                    E_s = self.galois_conjugate(s).change_ring(Kd)
                    psi = Ed.isogeny(x - l)
                    E_t = psi.codomain()
                    # Making sure the isomorphism is defined over Kd,
                    # extending Kd if necessary
                    c4s, c6s = E_s.c_invariants()
                    c4t, c6t = E_t.c_invariants()
                    if j_t == Integer(0) :
                        m, um = Integer(6) , c6t/c6s
                    elif j_t == Integer(1728) :
                        m, um = Integer(4) , c4t/c4s
                    else:
                        m, um = Integer(2) , (c6t*c4s)/(c6s*c4t)
                    f_iso = (x**m - um).factor()[Integer(0) ][Integer(0) ]
                    if f_iso.degree() > Integer(1) :
                        K_iso = Kd.extension(f_iso, names=('lu',)); (lu,) = K_iso._first_ngens(1)
                        Ed_iso = Ed.change_ring(K_iso)
                        E_s = E_s.change_ring(K_iso)
                        psi = Ed_iso.isogeny((x - l).change_ring(K_iso))
                        E_t = psi.codomain()
                    psi.set_post_isomorphism(E_t.isomorphism_to(E_s))
                    self._add_isogeny(s, psi)

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
        R = self.definition_field()
        Rx = PolynomialRing(R, names=('x',)); (x,) = Rx._first_ngens(1)
        Fx = Rx.fraction_field()
        phi_x = self._phi_x[sigma]
        num = phi_x.numerator()
        den = phi_x.denominator()
        _, iota = write_as_extension(self._to_Kphi, give_map=True)
        result = (Fx(sum(iota(num.coefficient(e)).list()[Integer(0) ] * x**e[Integer(0) ]
                         for e in num.exponents())) /
                  Fx(sum(iota(den.coefficient(e)).list()[Integer(0) ] * x**e[Integer(0) ]
                         for e in den.exponents())))
        if change_ring == None:
            return result
        Rx2 = Rx.change_ring(change_ring.codomain())
        Fx2 = Rx2.fraction_field()
        change_ring = Fx.Hom(Fx2)(Rx.Hom(Rx2)(change_ring))
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
        Rx = PolynomialRing(L, names=["x"]).fraction_field()
        x = Rx.gens()[Integer(0) ]
        return {s : (Rx(F(s).numerator().change_ring(iota) /
                        F(s).denominator().change_ring(iota)),
                     from_Kphi(l(s)))
                for s in G}

    def base_extend(self, R):
        result = EllipticCurve_number_field.base_extend(self, R)
        if isinstance(result, EllipticCurve_number_field):
            K = self.definition_field()
            L = result.base_ring()
            if (hasattr(R, "codomain") and
                R.domain() == K and
                R.codomain() == L):
                K_to_L = R
            else:
                K_to_L = K.embeddings(L)
                if len(K_to_L) > Integer(0) :
                    K_to_L = K_to_L[Integer(0) ]
                else:
                    K_to_L = None
            if K_to_L != None:
                return Qcurve(result, isogenies=self._isogeny_data(K_to_L))
        return result

    def _minimal_a_invariants(self, from_min):
        r"""Give the a_invariants over the minimal definition_field

        INPUT:

        - ``from_min`` -- The map from the minimal definition field to
          the definition field.

        """
        _, iota = write_as_extension(from_min, give_map=True)
        return [iota(a).list()[Integer(0) ] for a in self.a_invariants()]
        
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
            R = Fsnum.parent().change_ring(Kmin)
            Fsnum = sum(iota(Fsnum.monomial_coefficient(m)).list()[Integer(0) ] *
                        R(m) for m in Fsnum.monomials())
            Fsden = sum(iota(Fsden.monomial_coefficient(m)).list()[Integer(0) ] *
                        R(m) for m in Fsden.monomials())
            Fsnum = Fsnum.change_ring(min_map)
            Fsden = Fsden.change_ring(min_map)
            S = R.change_ring(min_map.codomain()).fraction_field()
            result[s] = S(Fsnum / Fsden)
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
        return Qcurve(ainvs, isogenies=isogenies)
    
    def twist(self, gamma):
        r"""Give the twist of this Q-curve by a given element gamma.

        If this Q-curve was given by .. MATH::

            E : y^2 = x^3 + a_2 x^2 + a_4 x + a_6

        the twisted Q-curve is given by .. MATH::
        
            E : y^2 = x^3 + \gamma a_2 x^2 + \gamma^2 a_4 x
                      + \gamma^3 a_6

        INPUT:

        - ``gamma`` -- An element of a number field.

        OUTPUT:
        
        A Q-curve which is the twist of this Q-curve by gamma. This
        curve will be defined over the smallest possible field over
        which it is completely defined.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2]); E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 5 with t = 2.236067977499790?
            sage: E2 = E.twist(t); E2
            Q-curve defined by y^2 = x^3 + (-6*lu0)*x^2 + (-45*lu0+90)*x over Number Field in lu0 with defining polynomial x^2 - 20 with lu0 = -1/7*lu^3 + 13/7*lu
            sage: E2.complete_definition_field()
            Number Field in agamma00 with defining polynomial x^4 - 16*x^3 - 8*x^2 + 576*x - 1264

        """
        ainvs, isogenies = self._twist(gamma)
        return Qcurve(ainvs, isogenies=isogenies).minimize_fields()

    @cached_method
    def conductor_restriction_of_scalars(self):
        r"""Give the conductor of the restriction of scalars of this Q-curve.

        OUTPUT:

        The conductor of the restriction of scalars of this curve over
        the decomposition field.

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.conductor_restriction_of_scalars()
            5566277615616

        .. SEE_ALSO::

            :meth:`decomposition_field`

        """
        K0 = self.definition_field()
        K = self.decomposition_field()
        if K0 != K:
            return self._over_Kdec().conductor_restriction_of_scalars()
        else:
            # Proposition 1 of Milne, On the arithmetic of Abelian varieties
            return self.conductor().absolute_norm() * K.discriminant()**Integer(2) 
    
    def trace_of_frobenius(self, prime, power=Integer(1) , splitting_map=Integer(0) ):
        r"""Compute the trace of a Frobenius element under the Galois
        representation associated with this curve.

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

        OUTPUT:

        An algebraic number that, when cast to the correct field, is
        the trace of a Frobenius element at `prime` to the power
        `power` under the $\lambda$-adic or mod $\lambda$ Galois
        representation associated to the splitting map given by
        `splitting_map` for $\lambda$ not dividing `prime`.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-2)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: E.trace_of_frobenius(5)
            a
            sage: E.trace_of_frobenius(11)
            2

        Works in case the reduction is multiplicative at the given prime.

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-1)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t*79/3), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: K = E.definition_field()
            sage: P = K.prime_above(5)
            sage: E.has_multiplicative_reduction(P)
            True
            sage: E.trace_of_frobenius(5)
            -6

        The trace depends on the splitting map used::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: E.trace_of_frobenius(7, splitting_map=0)
            -1/3*zeta4a^3 - 1/3*zeta4a
            sage: E.trace_of_frobenius(7, splitting_map=1)
            1/3*zeta4a^3 + 1/3*zeta4a
            sage: E.trace_of_frobenius(7, splitting_map=2)
            -1/3*zeta4a^3 - 1/3*zeta4a
            sage: E.trace_of_frobenius(7, splitting_map=3)
            1/3*zeta4a^3 + 1/3*zeta4a

        A sufficiently high power of the trace corresponds to the
        trace of Frobenius from the normal galois representation of
        the curve::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: L = E.definition_field()
            sage: L == E.decomposition_field()
            True
            sage: P = L.prime_above(17)
            sage: F = P.residue_field()
            sage: f = F.degree()
            sage: 1 + len(F) - E.reduction(P).count_points() # Trace of Frob_P
            322
            sage: E.trace_of_frobenius(17, power=f)
            322

        TESTS::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: K = E.definition_field()
            sage: def check(p):
            ....:     P = K.prime_above(p)
            ....:     F = P.residue_field()
            ....:     f = F.degree()
            ....:     En = E.reduction(P).count_points()
            ....:     a1 = 1 + len(F) - En
            ....:     a2 = E.trace_of_frobenius(p, power=f)
            ....:     return a1 == a2
            ....: 
            sage: all(check(p) for p in prime_range(7, 100))
            True

        """
        if power > Integer(1) :
            T = self.trace_of_frobenius(prime,
                                        splitting_map=splitting_map)
            D = self.determinant_of_frobenius(prime,
                                              splitting_map=splitting_map)
            F = self._trace_power_formula(power)
            iota = self._splitting_field_map(splitting_map)
            return F(T, iota(D))

        if not self.does_decompose():
            raise ValueError("This Q-curve must decompose for " +
                             "the Galois representations over the " +
                             "rationals to extend the Galois " +
                             "representations over the decomposition " +
                             "field.")
        K = self.decomposition_field()
        if self.definition_field() != K:
            return self._over_Kdec().trace_of_frobenius(prime, power=power,
                                                        splitting_map=splitting_map)
        if prime.divides(K.absolute_discriminant()):
            print("Warning: The decomposition field is ramified " +
                  "at " + str(prime) + " and the Galois " +
                  "representations over the rationals are therefore " +
                  "probably not unramified at " + str(prime))
        beta = self.splitting_map(splitting_map)
        for P in K.primes_above(prime):
            if self.has_good_reduction(P):
                result = self._trace_of_frob_good(P, beta)
                if result != None:
                    return result
            if self.has_multiplicative_reduction(P):
                return self._trace_of_frob_mult(P, beta)
        raise ValueError("This curve does not have the right " +
                         "type of reduction at " + str(prime) +
                         " or the reduction of the isogeny " +
                         "is not separable.")

    def _trace_of_frob_good(self, P, beta):
        r"""Implementation of meth:`trace_of_frobenius` in case of good reduction"""
        # Setting up the isogeny between the right minimal models
        K = self.definition_field()
        Emin = self.local_data(P).minimal_model()
        phi_min = Emin.isomorphism_to(self)
        phi_min = _rational_maps_of_isomorphism(phi_min)
        G = K.galois_group()
        Frob = G.artin_symbol(P)
        sE = self.galois_conjugate(Frob)
        sP = Frob(P)
        sEmin = sE.local_data(sP).minimal_model()
        sphi_min = sE.isomorphism_to(sEmin)
        sphi_min = _rational_maps_of_isomorphism(sphi_min)
        phi = self._phi_x[Frob], self._phi_y[Frob]
        phi = phi[Integer(0) ](phi_min), phi[Integer(1) ](phi_min)
        phi = sphi_min[Integer(0) ](phi), sphi_min[Integer(1) ](phi)

        # Checking whether the reduction of the isogeny is separable
        R = P.residue_field()
        Rx = R['x']; (x,) = Rx._first_ngens(1)
        Rx = Rx.fraction_field()
        phi0 = (phi[Integer(0) ].numerator().change_ring(R) /
                phi[Integer(0) ].denominator().change_ring(R))
        F = Rx(phi0(x, Integer(0) ))
        if F.derivative(x) == Integer(0) :
            return None # This will cause the function above to try
                        # other primes or fail

        # Defining variables needed for both p = 2 and p != 2
        phi1 = (phi[Integer(1) ].numerator().change_ring(R) /
                phi[Integer(1) ].denominator().change_ring(R))
        H = phi1(x, Integer(0) )
        G = phi1(x, Integer(1) ) - H
        Ered = Emin.reduction(P)
        p = P.smallest_integer()
        c1 = (F - x**p).numerator()

        # Computing number of points in the set
        # {P : phi P = Frob P}
        if p == Integer(2) :
            g = Ered.a1()*x + Ered.a3()
            h = (x**Integer(3)  + Ered.a2()*x**Integer(2)  + Ered.a4()*x +
                 Ered.a6())
            c3 = (g*G*h + g*G*H + G**Integer(2) *h + g**Integer(2) *H + h**Integer(2)  +
                  H**Integer(2) ).numerator()
            c4 = (g - G).numerator()
            gc13 = gcd(c1, c3)
            gc134 = gcd(gc13, c4)
            gc134g = gcd(gc134, g)
            num = (Integer(1)  + gc13.radical().degree() +
                   gc134.radical().degree() -
                   gc134g.radical().degree())
        else:
            R = (Integer(4) *x**Integer(3)  + Ered.b2()*x**Integer(2)  + Integer(2) *Ered.b4()*x +
                 Ered.b6())
            sEred = sEmin.reduction(sP)
            l = _scalar_of_rational_maps(phi0, phi1, Ered, sEred)
            c2 = (l * R**((p + Integer(1) )/Integer(2) ) -
                  F.derivative(x) * R).numerator()
            num = (Integer(1)  + Integer(2)  * gcd(c1, c2).radical().degree() -
                   gcd(c1, R).radical().degree())

        # Computing a_p(E)
        apE = F.numerator().degree() + p - num

        # The final result
        return beta(Frob)**(-Integer(1) ) * apE

    def _trace_of_frob_mult(self, P, beta):
        r"""Implementation of :meth:`trace_of_frobenius` in case of multiplicative reduction"""
        gamma = - self.c4() / self.c6()
        K = self.definition_field()
        L, iota, sqrtgamma = field_with_root(K, gamma,
                                             names='sqrtgamma',
                                             give_embedding=True)
        if not L.is_galois():
            L, iota2 = L.galois_closure(names='sqrtgamma', map=True)
            iota = _concat_maps(iota, iota2)
            sqrtgamma = iota2(sqrtgamma)
        Q = L.prime_above(iota(P))
        G = L.galois_group()
        FrobQ = G.artin_symbol(Q)
        FrobP = galois_field_restrict(FrobQ, K, embedding=iota)
        if self.degree_map(FrobP) != Integer(1) :
            return None # The degree of the isogeny should be square
                        # free, hence 1
        phix = self._phi_x[FrobP].numerator()
        phiy = self._phi_y[FrobP].numerator()
        x, y = phiy.parent().gens()
        u2 = phix.monomial_coefficient(x)
        u3 = phiy.monomial_coefficient(y)
        u = u3 / u2
        apE = QQ(iota(u)**(-Integer(1) ) * sqrtgamma * (FrobQ(sqrtgamma))**(-Integer(1) ))
        p = P.smallest_integer()
        return beta(FrobP)**(-Integer(1) ) * apE * QQ(Integer(1)  + p)
        
    def determinant_of_frobenius(self, prime, power=Integer(1) , splitting_map=Integer(0) ):
        r"""Compute the determinant of a Frobenius element acting on this
        curve.

        Given that this Q-curve decomposes over its decomposition
        field, one can associate to each splitting map $\beta$ a
        $\QQ$-simple abelian variety $A_\beta$ of $GL_2$-type which
        has this Q-curve as a 1-dimensional quotient. For each finite
        prime $\lambda$ in the image field of $\beta$, this $A_\beta$
        defines a 2-dimensional $\lambda$-adic Galois representation
        of the absolute Galois group of $\QQ$ that extends the
        $l$-adic Galois representation of this curve, where $l$ is the
        prime number below $\lambda$.

        Let $p$ be a prime number distinct from $l$ such that $p$ does
        not divide the conductor of the splitting character
        corresponding to $\beta$. In that case the determinant of the
        $\lambda$-adic Galois representation is
        unramified. Furthermore the determinant of the Frobenius
        element at $p$ can be explicitly computed and is an algebraic
        number that only depends on the splitting character
        corresponding to $\beta$ and the prime number $p$.

        This function computes the algebraic number that is the
        determinant of (a power of) a Frobenius element. It will first
        check that the give prime number `prime` satisfies the
        mentioned condition and raise a ValueError if this is not the
        case.

        INPUT:

        - ``prime`` -- A prime number such that `prime` does not
          divide the conductor of the splitting character associated
          to the splitting map corresponding to `splitting_map`.

        - ``power`` -- A strictly positive integer (default: 1). If
          set to a value higher than 1 will compute the trace of the
          Frobenius element to the given power instead of the
          Frobenius element itself.

        - ``splitting_map`` -- A non-negative integer smaller than the
          number of splitting maps associated to this Q-curve
          (default: 0). This indicates the splitting map for which the
          determinant of frobenius should be computed.

        OUTPUT:

        An algebraic number that is the determinant of the Frobenius
        element at `prime` to the power `power` under any
        $\lambda$-adic Galois representation associated to the
        splitting map given by `splitting_map` for $\lambda$ not
        dividing `prime`.

        EXAMPLES::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: E.determinant_of_frobenius(17)
            -17*zeta4

        The determinant will depend on the chosen splitting map::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: E.determinant_of_frobenius(7, splitting_map=0)
            7*zeta4
            sage: E.determinant_of_frobenius(7, splitting_map=4)
            -7*zeta4

        The determinant of a power of frobenius corresponds to the
        determinant of the standard galois representation associated
        with an elliptic curve::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-1)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E = E.decomposable_twist()
            sage: L = E.definition_field()
            sage: L == E.decomposition_field()
            True
            sage: P = L.prime_above(11)
            sage: F = P.residue_field()
            sage: len(F) # determinant of Frob_P
            121
            sage: E.determinant_of_frobenius(11, power=F.degree())
            121

        """
        if power > Integer(1) :
            D = self.determinant_of_frobenius(prime,
                                              splitting_map=splitting_map)
            return D**power
        else:
            eps = self.splitting_character(splitting_map)
            if prime.divides(eps.conductor()):
                raise ValueError("The given prime number " +
                                 str(prime) +
                                 " divides the conductor of the " +
                                 "splitting character corresponding "+
                                 "to the given splitting map.")
            return eps(prime)**(-Integer(1) ) * prime
            
    @cached_method
    def newform(self, algorithm='sage', verify=Integer(0) , warning=Integer(10) **Integer(4) ,
                path=None):
        r"""Give a newform associated to this Q-curve.

        Each non-CM Q-curve is the quotient of a $\Q$-simple variety
        of GL_2-type, which in turn is isogeneous to an abelian
        variety associated to a newform. We will call such a newform a
        newform associated to a Q-curve.

        ALGORITHM:

        To obtain the newforms we take the restriction of scalars of
        this Q-curve over its decomposition field. If all is well,
        this is isogenous to a product of $\Q$-simple, non
        $\Q$-isogenous abelian varieties of GL_2-type. The newforms
        attached to these varieties are twists of one another,
        allowing the computation of their possible levels by using
        formulas linking the levels of twisted newforms, the conductor
        of the restriction of scalars and the conductor of the abelian
        varieties.
        
        Next for each possible combination of levels, the space of
        newforms of the lowest level is computed with the appropriate
        character. For these newforms the traces of Frobenius of this
        Q-curve and the traces of Frobenius of the newforms are
        compared until a single newform remains. Some Q-curves have
        multiple factors of the same level and character, so one is
        chosen.

        INPUT:
        
        - ``algorithm`` -- A string that determines which program
          should be used to compute the spaces of newforms. Allowed
          values are: 'sage' (default) to use Sage, 'magma' to use
          MAGMA, or 'file' to load the newforms from a file.

        - ``verify`` -- A non-negative integer determining what the
          biggest prime is for which the result should be verified
          using the traces of the corresponding Frobenius elements.

        - ``warning`` -- A non-negative integer (default 10^4). If the
          primes at which traces of frobenius are compared ever exceed
          this value, a warning will be printed to the standard output
          and this value will be doubled. If it is smaller than the
          value of `verify` it will be set to `verify`.

        - ``path`` -- A string or None (default: None). Only used in
          case `algorithm` is set to 'file', in which case this should
          be the (relative) path to the file from which the newforms
          should be loaded as a string.

        OUTPUT:

        A tuple consisting of

        - a wrapped newform such that this curve is a quotient of the
          abelian variety associated to that newform, i.e. a newform
          associated to this Q-curve by modularity

        - A list of Dirichlet characters that twist this newform into
          other newforms, such that the restriction of scalars of this
          Q-curve over the decomposition field is isogenous (over
          $\Q$) to the product of the abelian varieties associated to
          the twists of the first newform by the given characters.

        .. NOTE::
        
        The newforms associated to a Q-curve are in no way unique. Not
        only are all their galois conjugates also associated, also
        twists of the newform are associated to the same curve. On the
        other hand, all curves isogenous to this Q-curve are
        associated to the same collection of newforms.

        .. SEEALSO::

            :meth:`decomposition_field`,
            :meth:`does_decompose`,
            :meth:`splitting_character`,
            :meth:`twist_character`,
            :meth:`conductor_restriction_of_scalars`

        EXAMPLE::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E2 = E.decomposable_twist()
            sage: E2.newform() # long time (22 minutes)
            (q - 1/2*a3*q^3 + (-1/2*a3 - 1)*q^5 + O(q^6),
             [Dirichlet character modulo 1 of conductor 1])

        Or another example with two factors::

            sage: from modular_method.elliptic_curves.Qcurves import Qcurve
            sage: K.<t> = QuadraticField(-3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E2 = E.decomposable_twist()
            sage: E2.newform() # long time (12 seconds)
            (q + a0*q^3 + (-a0^2 + 1)*q^5 + O(q^6),
             [Dirichlet character modulo 8 of conductor 8 mapping 7 |--> 1, 5 |--> -1,
              Dirichlet character modulo 1 of conductor 1])

        """
        if warning < verify:
            warning = verify
        if not self.does_decompose():
            raise ValueError("Can not compute newform if the restriction of " +
                             "scalars does not decompose.")
        
        levels = self.newform_levels()
        twists_base = self.twist_character('conjugacy')
        # Find a common base for all twists
        M = lcm(chi.modulus() for chi in twists_base)
        twists_base = [chi.extend(M)**(-Integer(1) ) for chi in twists_base]
        L = QQ
        for i in range(len(twists_base)):
            Li = twists_base[i].base_ring()
            if L != Li:
                L, old_to_new, i_to_new = composite_field(L, Li,
                                                          give_maps=True)
                twists_base[i] = twists_base[i].change_ring(i_to_new)
                for j in range(i):
                    twists_base[j] = twists_base[j].change_ring(old_to_new)
        # Characters of the newforms
        eps_ls = [(eps**(-Integer(1) )).primitive_character()
                  for eps in self.splitting_character('conjugacy')]
        Lbeta_ls = self.splitting_image_field('conjugacy') # coefficient fields

        if (algorithm != 'magma' and algorithm != 'sage' and
            algorithm != 'file'):
            raise ValueError("%s is not a valid algorithm to use."%algorithm)

        candidates = []
        max_level = lcm(lcm(tmp) for tmp in levels)
        # Keeps track of the lcm of all N considered
        # the primes in these will be excluded in checking
        # the traces of Frobenius
        done_cases = {}
        done_cases2 = []
        for k in range(len(levels)):
            # Newform with smallest level:
            i_min, N = min(enumerate(levels[k]), key=(lambda x: x[Integer(1) ]))
            eps = eps_ls[i_min]
            if (N, eps, i_min) not in done_cases2:
                chi = twists_base[i_min]
                # Twists relative to i_min
                twists = [(chi_j * chi**(-Integer(1) )).primitive_character().minimize_base_ring()
                          for chi_j in twists_base]
                Lbeta = Lbeta_ls[i_min]

                # Computing the newforms, but only if not done already
                if (N, eps) not in done_cases:
                    done_cases[(N, eps)] = get_newforms(N, character=eps,
                                                        algorithm=algorithm,
                                                        path=path)
                nfs = done_cases[(N, eps)]
                for nf in nfs:
                    Kf = nf.coefficient_field()
                    if Kf.absolute_degree() == Lbeta.absolute_degree():
                        for iota in Kf.embeddings(Lbeta):
                            candidates.append((nf, twists, i_min, iota))
                # Making sure we don't do this twice            
                done_cases2.append((N, eps, i_min))

        p = Integer(1) 
        while len(candidates) > Integer(1)  or p < verify:
            if p > warning:
                print("Warning: Checked prime exceeds warning value " +
                       str(warning))
                warning = Integer(2)  * warning
            while p.divides(max_level):
                p = next_prime(p)
            removed = []
            for nf, twists, i_min, iota in candidates:
                apf = iota(nf.trace_of_frobenius(p))
                splitting_map = self._conjugacy_classes()[i_min][Integer(0) ]
                apE = self.trace_of_frobenius(p, splitting_map=splitting_map)
                if apf != apE:
                    removed.append((nf, twists, i_min, iota))
            for candidate in removed:
                candidates.remove(candidate)
            p = next_prime(p)

        if len(candidates) < Integer(1) :
            raise ValueError("No newform seems to correspond with " +
                             "this curve. A bug?")
        return candidates[Integer(0) ][Integer(0) ], candidates[Integer(0) ][Integer(1) ]

    def _repr_(self):
        """
        String representation of a Q-curve.

        REMARK:

        This is a direct copy from the code included
        in EllipticCurve_number_field

        """
        b = self.ainvs()
        a = [z._coeff_repr() for z in b]
        s = "Q-curve defined by "
        s += "y^2 "
        if a[Integer(0) ] == "-1":
            s += "- x*y "
        elif a[Integer(0) ] == '1':
            s += "+ x*y "
        elif b[Integer(0) ]:
            s += "+ %s*x*y "%a[Integer(0) ]
        if a[Integer(2) ] == "-1":
            s += "- y "
        elif a[Integer(2) ] == '1':
            s += "+ y "
        elif b[Integer(2) ]:
            s += "+ %s*y "%a[Integer(2) ]
        s += "= x^3 "
        if a[Integer(1) ] == "-1":
            s += "- x^2 "
        elif a[Integer(1) ] == '1':
            s += "+ x^2 "
        elif b[Integer(1) ]:
            s += "+ %s*x^2 "%a[Integer(1) ]
        if a[Integer(3) ] == "-1":
            s += "- x "
        elif a[Integer(3) ] == '1':
            s += "+ x "
        elif b[Integer(3) ]:
            s += "+ %s*x "%a[Integer(3) ]
        if a[Integer(4) ] == '-1':
            s += "- 1 "
        elif a[Integer(4) ] == '1':
            s += "+ 1 "
        elif b[Integer(4) ]:
            s += "+ %s "%a[Integer(4) ]
        s = s.replace("+ -","- ")
        s += "over %s"%self.definition_field()
        return s

