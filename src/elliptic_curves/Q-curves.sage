r"""A class for working with $\Q$-curves

A $\Q$-curve is an elliptic curve defined over some number field that
is isogenous to all its galois conjugates. This file contains the
class Qcurve that represents $\Q$-curves that do not have complex
multiplication.

EXAMPLES::

    sage: K.<t> = QuadraticField(-2)
    sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : (-t, 2)})
    Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2

Q-curves without CM are modular and are linked to classical newforms
that can be computed::

    sage: K.<t> = QuadraticField(3)
    sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
    sage: E2 = E.decomposable_twist()
    sage: E2.newform() # long
    (q + (-a + 1)*q^3 + a*q^5 + 3*a*q^7 + (-2*a - 1)*q^9 + 4*q^11 + O(q^12),
     [Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1])

AUTHORS:

- Joey van Langen (2019-03-01): initial version

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
from sage.schemes.elliptic_curves.ell_number_field import EllipticCurve_number_field

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

        sage: E = EllipticCurve([0,0,0,1,0])
        sage: P = E.torsion_points()[0]
        sage: phi = E.isogeny(P)
        sage: _scalar_of_isogeny(phi)
        1

    Note that the scalar of an isogeny and its dual multiply to the
    degree of the isogeny::

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
    x,y = Fy.parent().gens()
    f1 = phi.domain().defining_polynomial()(x, y, 1)
    f2 = phi.codomain().defining_polynomial()(x, y, 1)
    R = Fy.parent().base().base()
    u = R(((Fx.derivative(x) * f1.derivative(y)) / f2.derivative(y)(Fx, Fy)).numerator())
    f = u.minpoly()
    return NumberField(f, names='l').gen()

def _scalar_of_isomorphism(E1, E2):
    """Return the scalar associated to an isomorphism.
    
    For an isomorphism of the form .. MATH:

        (x,y) \mapsto (F(x), y F'(x) / \lambda )

    returns $ 1 / u $.

    If both elliptic curves are defined over a subfield of the complex
    numbers, this scalar is the same as the scalar in the map ..MATH

       z \mapsto \lambda z

    on the complex numbers that induces this isogeny on the
    corresponding quotients.

    INPUT:

    - ``E1`` -- An elliptic curve given by a Weierstrass equation over
      a field of characteristic not equal to 3 and given by a
      Weierstrass equation of the form .. MATH:
    
        Y^2 = X^3 + a2*X^2 + a4*X + a6

    - ``E2`` -- An elliptic curve defined over the same basefield as
      E1, which is isomorphic to E1 and the isomorphism can be given
      in a form as above. The isomorphism does not have to be defined
      over the base field.

    OUTPUT:

    The unique number $\lambda$ such that the isomorphism from E1 to
    E2 is of the form .. MATH:

        (x,y) \mapsto (F(x), y F'(x) / \lambda)
    
    or an error if such a $\lambda$ does not exist.

    """
    if (E1.a1() != 0 or E1.a3() != 0):
        raise ValueError("The curve %s is not of the correct form."%(E1,))
    if (E2.a1() != 0 or E2.a3() != 0):
        raise ValueError("The curve %s is not of the correct form."%(E2,))
    E1 = E1.rst_transform(-E1.a2()/3,0,0);
    E2 = E2.rst_transform(-E2.a2()/3,0,0);
    ainv1 = list(E1.a_invariants())
    ainv1.insert(0,0); ainv1.insert(5,0);
    ainv2 = list(E2.a_invariants())
    ainv2.insert(0,0); ainv2.insert(5,0);
    k = []
    u = []
    for i in range(len(ainv1)):
        if ainv1[i] != 0 and ainv2[i] != 0:
            k.append(i)
            u.append(ainv2[i]/ainv1[i])
        else:
            u.append(0)
    n = gcd(k)
    if n in k:
        un = u[n]
    elif n + min(k) in k:
        un = u[n + min(k)]/(u[min(k)])
    else:
        raise Exception("Could not compute the sought scalar for some reason.")
    if hasattr(un, "nth_root"):
        try:
            return 1 / un.nth_root(n)
        except ValueError:
            pass
    K = un.parent()
    R.<x> = K[]
    L.<u> = K.extension(x^n - un)
    return 1 / L.absolute_field(names='u')(u)

class Qcurve(EllipticCurve_number_field):
    r"""A Q-curve over some number field

    A Q-curve is an elliptic curve defined over some number field,
    such that all its galois conjugates are isogeneous to the curve
    itself.

    In this class a Q-curve is represented as an elliptic curve E
    defined over a galois number field K, together with for each
    element s of the galois group of K a tuple containing the scalar
    and degree associated to the isogeny from s(E) to E.

    .. NOTE::

    This class is intended for Q-curves without complex
    multiplication. Although Q-curves with complex multiplication
    might work, the theory behind many of the methods in this class
    was only intended for Q-curves without complex multiplication.

    EXAMPLE::

        sage: K.<t> = QuadraticField(-2)
        sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : (-t, 2)})
        Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2

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

        Third it will fill in data about the isogenies from a galois
        conjugate to itself using the provided data given.

        Next it will determine the remaining isogenies using the
        guessed degrees and combining previously found isogenies.
        
        If in this way it can not find isogenies from each galois
        conjugate to itself, the initialization will produce an error
        and the resulting object should not be used.

        INPUT:
        
         - ``curve`` -- An elliptic curve over some Galois number
           field or any input that would create such a curve when
           passed to the constructor EllipticCurve. It should be
           defined by a Weierstrass equation of the form .. MATH:

              Y^2 + a_1 X Y + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6

         - ``isogenies`` -- A dictionary (default: {}) with as keys
           elements of the galois group of the base field of the
           Q-curve and as values the corresponding isogeny from the
           galois conjugate of this Q-curve to itself. Such isogenies
           must be defined over the field over which the curve is
           defined or a direct extension thereof. The isogeny can be
           given as either an isogeny as a Sage object; a tuple of a
           rational function in $x$ and an algebraic number, that are
           respectively the $x$-coordinate map of the isogeny and the
           induced scalar multiplication on the differentials; or a
           tuple of three rational functions $F$, $G$, $H$ in $x$ such
           that the isogeny is $(x, y) \mapsto (F(x), G(x) y + H(x))$
           outside points mapping to infinity.

        - ``guessed_degrees`` -- A list (default: []) of strictly
           positive integers indicating possible degrees of isogenies
           from galois conjugates of this curve to itself.

        - ``verbose`` -- A boolean value or an integer (default:
           False). When set to True or any value larger then zero will
           print comments to stdout about the computations being done
           whilst busy. If set to False or 0 will not print such
           comments. If set to any negative value will also prevent
           the printing of any warnings. A higher value will cause
           more messages to be printed.

        EXAMPLES:

        One can create a Qcurve using an explicit isogeny::

            sage: K.<t> = QuadraticField(-2)
            sage: G.<s> = K.galois_group()
            sage: E = EllipticCurve([0, 12, 0, 18*(t + 1), 0])
            sage: Es = EllipticCurve([s(a) for a in E.a_invariants()])
            sage: phi = Es.isogeny(Es([0,0]), codomain=E)
            sage: Qcurve(E, isogenies={s : phi})
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2

        However, it is sufficient to simply provide the scalar and the
        degree of the isogeny::

            sage: K.<t> = QuadraticField(-2)
            sage: Qcurve([0, 12, 0, 18*(t + 1), 0], isogenies={s : (-t, 2)})
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 + 2
        
        It is also possible to make the code 'guess' the isogenies by
        giving a suggestion for their degree::
        
            sage: K.<t> = QuadraticField(3)
            sage: Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 3

        Note that giving no data with regards to the isogenies will
        result in an error::

            sage: K.<t> = QuadraticField(5)
            sage: Qcurve([0, 12, 0, 18*(t + 1), 0])
            Traceback (most recent call last):
            ...
            ValueError: There is not sufficient isogeny information to make [0, 12, 0, 18*t + 18, 0] a Q-curve

        """
        self._init_curve(curve)
        self._init_isogenies()
        for sigma, phi in isogenies.iteritems():
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
        if not all(d(s) * d(t) / d(s*t) == c(s,t)^2 for s in G for t in G):
            raise ValueError("The given isogenies are not valid.")

    def definition_field(self):
        r"""Give the field over which this Q-curve is defined.

        .. NOTE::

        This number field is always galois.

        OUTPUT:

        The number field over which this Q-curve is defined.

        EXAMPLE::

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.definition_field()
            Number Field in t with defining polynomial x^2 - 3

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

    def _galois_cache_key(self, sigma):
        r"""Give a cache key for an element of a galois group"""
        return str(sigma), sigma.parent().number_field()
        
    @cached_method(key=_galois_cache_key)
    def galois_conjugate(self, sigma):
        r"""Give the galois conjugate of this curve.

        INPUT:

        - ``sigma`` -- A galois homomorphism of some number field

        OUTPUT:
        
        The galois conjugate of this curve by the galois homomorphism
        which extends to a common galois homomorphism over the
        algebraic closure of Q as sigma. This will be an elliptic
        curve and not a Q-curve

        EXAMPLE::

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 3
            sage: sigma = K.galois_group().gens()[0]
            sage: E.galois_conjugate(sigma)
            Elliptic Curve defined by y^2 = x^3 + 12*x^2 + (-18*t+18)*x over Number Field in t with defining polynomial x^2 - 3

        """
        sigma = galois_field_change(sigma, self.definition_field())
        return EllipticCurve(self.a_invariants()).change_ring(sigma.as_hom())

    # Isogeny related stuff
    def _init_isogenies(self):
        r"""Initialize the isogeny data.

        """
        # The x & y coordinate maps of isogenies
        self._phi_x = dict()
        self._phi_y = dict()
        # Map from the base field of the elliptic curve:
        self._Kphi = self.definition_field()
        # Initialize the trivial isogeny that is there:
        R = self.definition_field()
        e = R.galois_group().identity()
        Rxy = PolynomialRing(R, names=["x", "y"]).fraction_field()
        x, y = Rxy.gens()
        self._phi_x[e] = x
        self._phi_y[e] = y

    def _add_isogeny(self, sigma, phi):
        r"""Add an isogeny to the stored isogeny data.

        INPUT:

        - ``sigma`` -- A galois homomorphism of the field over which
          this Q-curve is defined.

        - ``phi`` -- An isogeny from the galois conjugate of this
          curve by sigma to this curve itself, a tuple of the
          corresponding x-coordinate map and a scalar, or a tuple of
          the rational functions in x defining the isogeny.

        """
        if isinstance(phi, tuple):
            R = phi[0].base_ring()
            Rxy = PolynomialRing(R, names=['x','y']).fraction_field()
            x, y = Rxy.gens()
            self._phi_x[sigma] = phi[0](x)
            if len(phi) == 2:
                dom = self.galois_conjugate(sigma)
                codom = self
                self._phi_y[sigma] = ((self._phi_x[sigma].derivative(x)
                                       * (2*y + dom.a1()*x + dom.a3())
                                       / R(phi[1]))
                                      - (codom.a1()*self._phi_x[sigma] + codom.a3())) / 2
            elif len(phi) == 3: 
                self._phi_y[sigma] = phi[1](x)*y + phi[2](x)
        else:
            self._add_isogeny(sigma, (phi.x_rational_map(), _scalar_of_isogeny(phi)))

    def _has_isogeny(self, sigma):
        r"""Tell if an isogeny from a Galois conjugate is already registered.

        """
        return (sigma in self._phi_x and
                self._phi_x[sigma] != None and
                sigma in self._phi_y and
                self._phi_y[sigma] != None)

    def _update_isogeny(self, sigma, change):
        r"""Update the field of one specific isogeny"""
        self._phi_x[sigma] = self._phi_x[sigma].change_ring(change)
        self._phi_y[sigma] = self._phi_y[sigma].change_ring(change)

    def _update_isogeny_field(self):
        r"""Update the field over which all isogenies are defined

        """
        G = list(self.definition_field().galois_group())
        for i in range(len(G)):
            if self._has_isogeny(G[i]):
                Ri = self._phi_x[G[i]].base_ring()
                if Ri != self._to_Kphi.codomain():
                    _, old_to_new, Ri_to_new = composite_field(self._Kphi, Ri,
                                                               give_maps=True)
                    def_to_new = (old_to_new *
                                  self.definition_field().hom(self._Kphi))
                    self._Kphi, iota = write_as_extension(def_to_new,
                                                          give_map=True)
                    old_to_new = iota * old_to_new
                    Ri_to_new = iota * Ri_to_new
                    self._update_isogeny(G[i], Ri_to_new)
                    for j in range(i):
                        if _has_isogeny(self, G[j]):
                            self._update_isogeny(G[j], old_to_new)
        if not self._Kphi.is_galois():
            Kphi.<al> = self._Kphi.galois_closure()
            clos = self._Kphi.embeddings(Kphi)[0]
            def_to_clos = clos * self.definition_field().hom(self._Kphi)
            self._Kphi, old_to_new = write_as_extension(def_to_clos,
                                                        give_map=True)
            for s in G:
                if self._has_isogeny(G[i]):
                    self._update_isogeny(G[i], old_to_new)

    def _fill_isogenies(self):
        r"""Attempt to fill in missing isogenies by combining known ones.

        """
        self._update_isogeny_field()
        G = self.definition_field().galois_group()
        for s in G:
            for t in G:
                if (not self._has_isogeny(s*t) and
                    self._has_isogeny(s) and
                    self._has_isogeny(t)):
                    tL = galois_field_extend(t, self._Kphi)
                    tL_phi_s = (self._phi_x[s].change_ring(tL.as_hom()),
                                self._phi_y[s].change_ring(tL.as_hom()))
                    self._phi_x[s*t] = self._phi_x[t](tL_phi_s)
                    self._phi_y[s*t] = self._phi_y[t](tL_phi_s)
        return all(self._has_isogeny(s) for s in G)

    def _add_isogenies_of_degree(self, degree, verbose=False):
        r"""Attempt to find isogenies of a given degree.

        """
        G = self.definition_field().galois_group()
        fd = self.torsion_polynomial(degree)
        for g,e in fd.factor():
            Kd.<l> = self.definition_field().extension(g)        
            Ed = EllipticCurve_number_field.base_extend(self, Kd)
            S.<x> = Kd[]
            psi = Ed.isogeny(x - l)
            E_t = psi.codomain()
            j_t = E_t.j_invariant()
            for s in G:
                if Kd(self.galois_conjugate(s).j_invariant()) == j_t:
                    if verbose > 0:
                        print "Degree %s isogeny found for"%degree, s
                    E_s = self.galois_conjugate(s).change_ring(Kd)
                    psi.set_post_isomorphism(E_t.isomorphism_to(E_s))
                    self._add_isogeny(s, psi.dual())

    @cached_method(key=_galois_cache_key)
    def isogeny_scalar(self, sigma):
        r"""Return the scalar of the isogeny from the sigma conjugate of this
        curve to this curve.

        The scalar of an isogeny is the constant $\lambda$ such that
        the isogeny becomes $z \mapsto \lambda z$ on the complex
        numbers.

        INPUT:
        
        - ``sigma`` -- A galois homomorphism of a number field

        OUTPUT:

        The constant $\lambda$ such that the map $z \mapsto \lambda z$
        on the complex number induces an isogeny from a galois
        conjugate of this curve to itself. The galois conjugate is one
        obtained by conjugating with an extension of sigma.

        EXAMPLE::

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: G = K.galois_group()
            sage: E.isogeny_scalar(G[0])
            1
            sage: E.isogeny_scalar(G[1])
            1/10*tu0^3 + 3/10*tu0
            sage: E.isogeny_scalar(G[1]).parent()
            Number Field in tu0 with defining polynomial x^4 - 2*x^2 + 25

        """
        sigma = galois_field_change(sigma, self.definition_field())
        Fx, Fy = (self._phi_x[sigma], self._phi_y[sigma])
        x,y = Fy.parent().gens()
        f1 = self.galois_conjugate(sigma).defining_polynomial()(x, y, 1)
        f2 = self.defining_polynomial()(x, y, 1)
        R = Fy.parent().base_ring()
        return R(((Fx.derivative(x) * f1.derivative(y)) /
                  f2.derivative(y)(Fx, Fy)).numerator())

    @cached_method(key=_galois_cache_key)
    def isogeny_x_map(self, sigma):
        r"""Return the x-coordinate rational map of the isogeny from the sigma
        conjugate of this curve to this curve.

        The x-coordinate of the image of a point under an isogeny can
        be described as a rational function of the x-coordinate of the
        corresponding point in the domain.

        INPUT:
        
        - ``sigma`` -- A galois homomorphism of a number field

        OUTPUT:

        A rational function in $x$ that gives the $x$-coordinate of an
        image point of the isogeny as a rational function in the $x$
        coordinate of the origin point. The isogeny is the registered
        isogeny from the `sigma` galois conjugate of this curve to
        this curve.

        EXAMPLE::

        TODO

        """
        sigma = galois_field_change(sigma, self.definition_field())
        return self._phi_x[sigma]

    def complete_definition_field(self):
        r"""Give the field over which the Q-curve is completely defined.

        OUTPUT:

        A number field over which both this elliptic curve and all
        isogenies from its galois conjugates to itself are defined.

        EXAMPLE::

            sage: K.<t> = QuadraticField(-2)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.complete_definition_field()
            Number Field in t with defining polynomial x^2 + 2

        In general this field is bigger than the definition field::

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.definition_field()
            Number Field in t with defining polynomial x^2 - 3
            sage: E.complete_definition_field()
            Number Field in tu0 with defining polynomial x^4 - 2*x^2 + 25

        """
        return self._Kphi

    @cached_method(key=_galois_cache_key)
    def degree_map(self, sigma):
        r"""Give the degree of an isogeny from the sigma galois conjugate of
        this curve to this curve.

        INPUT:

        - ``sigma`` -- A galois homomorphism of a number field

        OUTPUT:
        
        The degree of an isogeny from a galois conjugate of this curve
        to this curve itself. The galois conjugate is one obtained by
        conjugating with an extension of the given galois homomorphism
        sigma.

        EXAMPLE::

            sage: K.<t> = QuadraticField(2)
            sage: G.<s> = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.degree_map(s)
            2
        
        The isomorphism from the curve to itself will always have a
        degree that is a square::

            sage: K.<t> = QuadraticField(3)
            sage: G = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.degree_map(G(1))
            1

        One can also use galois homomorphisms not necessarily defined
        over the definition field::

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
        return Fx.degree(x)

    @cached_method
    def degree_map_image(self):
        r"""Give the image of the degree map in $\Q^*/(\Q^*)^2$

        OUTPUT:
        
        A list of squarefree integers such that each value of the
        degree map is a square times such an integer and all integers
        in this list differ a square from a value of the degree map.

        EXAMPLE::

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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.degree_field()
            Number Field in t with defining polynomial x^2 - 3

        The degree field is always a subfield of the definition field,
        but can be strictly smaller::

            sage: K.<t> = CyclotomicField(3)
            sage: G.<s> = K.galois_group()
            sage: E = Qcurve([0, -6*t + 2, 0, t, 0], isogenies={s : (2, 4)})
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
        products = [1]
        if a1 != None and Kd(a1).is_square():
            ai.append(a1)
            products.append(a1)
        for tmp in Kd.subfields(degree=2):
            a = tmp[0].discriminant().squarefree_part()
            if a not in products:
                ai.append(a)
                products.extend([(a*b).squarefree_part() for b in products])
        di = [0]*len(ai)
        for sigma in Kd.galois_group():
            ls = [sigma(sqrt(Kd(a)))/sqrt(Kd(a)) for a in ai]
            if sum(ls) == len(ls)-2: # Precisely one entry == -1
                for i in range(len(ls)):
                    if ls[i] == -1:
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

            sage: K.<t> = QuadraticField(5)
            sage: G = K.galois_group()
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: matrix([[E.c(s, t) for t in G] for s in G])
            [ 1  1]
            [ 1 -2]

        Note that the value of this function always squares to the
        coboundary of the degree map::

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
        return QQ(l(sigma) * sigma(l(tau)) * l(sigma*tau)^(-1))

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
        if 2 not in result:
            result.insert(0,2)
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
            return 1
        else:
            return product([hilbert_symbol(ai,di,p)
                            for (ai,di) in self.xi_pm()])

    def _first_splitting_character(self):
        r"""Compute the first splitting character"""
        N = 1
        primes = [p for p in self._xi_pm_primes() if self.xi_pm_local(p) == -1]
        L = CyclotomicField(lcm(euler_phi(p) for p in primes))
        eps_ls = [DirichletGroup(1, base_ring=L)[0]]
        for p in primes:
            if p == 2:
                N *= 4
                eps_ls.append(DirichletGroup(4, base_ring=L).gen())
            else:
                N *= p
                eps_ls.append(DirichletGroup(p, base_ring=L).gen()^((p-1).odd_part()))
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
        if not self._is_cached('_eps') or 0 not in self._eps:
            self._eps = dict()
            self._eps[0] = [self._first_splitting_character()]
        if hasattr(i, "__iter__"):
            return tuple(self._splitting_character_data(ii, j) for ii in i)
        if i in ZZ:
            if i not in self._eps:
                eps0 = self._eps[0][0]
                chi = self.twist_character(i)
                N = lcm(eps0.conductor(),chi.conductor())
                L = composite_field(eps0.base_ring(), chi.base_ring())
                D = DirichletGroup(N, base_ring=L)
                epsi = D(eps0) * D(chi)^2
                epsi = epsi.primitive_character()
                epsi = epsi.minimize_base_ring()
                self._eps[i] = [epsi]
            if j >= 1 and len(self._eps[i]) < 2:
                self._eps[i].append(dirichlet_fixed_field(self._eps[i][0]))
            if j >= 2 and len(self._eps[i]) < 3:
                self._eps[i].append(dirichlet_to_galois(self._eps[i][0]))
            return self._eps[i][j]
        if i == 'all':
            return tuple(self._splitting_character_data(ii, j)
                         for ii in range(self.number_of_splitting_maps()))
        if i == 'conjugacy':
            return tuple(self._splitting_character_data(ii[0], j)
                         for ii in self._conjugacy_classes())
        raise Exception("Invalid index %s."%i)
    
    def splitting_character(self, index=0, galois=False):
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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(t + 1), 0], guessed_degrees=[2])
            sage: E.splitting_character()
            Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1

        There are as many splitting characters as the degree of the
        decomposition field::

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
            return self._splitting_character_data(index, 2)
        else:
            return self._splitting_character_data(index, 0)

    def splitting_character_field(self, index=0):
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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.splitting_character_field()
            Number Field in zeta0 with defining polynomial x^2 - 3

        In general it is distinct from the complete definition field::

            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: K1 = E.complete_definition_field(); K1
            Number Field in tu0 with defining polynomial x^4 - 6*x^2 + 49
            sage: K2 = E.splitting_character_field(); K2
            Number Field in zeta0 with defining polynomial x^4 - 5*x^2 + 5
            sage: K1.is_isomorphic(K2)
            False

        """
        return self._splitting_character_data(index, 1)

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
        if 2.divides(Keps.degree()):
            b = Keps.subfields(degree=2)[0][0].discriminant().squarefree_part()
        ai, di = self.dual_basis(a1=b)
        L = CyclotomicField(2*eps.order())
        for i in range(len(di)):
            Kdi = QuadraticField(di[i])
            if ai[i] == b:
                Lbig, L_to_Lbig, Kdi_to_Lbig = composite_field(L, Kdi,
                                                               give_maps=True)
                alpha = L_to_Lbig(L.gen()) * Kdi_to_Lbig(Kdi.gen())
                L = Lbig.subfield(alpha)[0]
            else:
                L = composite_field(L, Kdi)
        return L

    @cached_method
    def splitting_image_field(self, index=0):
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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12*t - 12, 0, 180 - 108*t, 0], guessed_degrees=[2])
            sage: K1 = E.splitting_image_field(); K1
            Number Field in zeta4a0 with defining polynomial x^2 + 2
            sage: beta = E.splitting_map()
            sage: G = E.splitting_field().galois_group()
            sage: [beta(s) for s in G]
            [1, zeta4a0]

        """
        eps = self.splitting_character(index)
        Keps = self.splitting_character_field(index)
        return self._splitting_image_field(eps, Keps)

    def _splitting_field(self, Keps):
        r"""Imlementation of :meth:`splitting_field`"""
        if isinstance(Keps, tuple):
            return tuple(self._splitting_field(Keps_i) for Keps_i in Keps)
        Kd = self.degree_field()
        return composite_field(Kd, Keps)

    @cached_method
    def splitting_field(self, index=0):
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

        OUTPUT:

        The splitting field corresponding to the splitting map of the
        given index. If the index was 'all' or 'conjugacy' will return
        a tuple of such fields corresponding to the corresponding
        tuple of indices. If the given index was a tuple will return a
        tuple of outputs on each entry of this tuple in the same
        order.

        EXAMPLE::

            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.splitting_field()
            Number Field in zeta0 with defining polynomial x^4 - 5*x^2 + 5
            sage: E.degree_field()
            Number Field in t with defining polynomial x^2 - 5
            sage: E.splitting_character_field()
            Number Field in zeta0 with defining polynomial x^4 - 5*x^2 + 5

        """
        Keps = self.splitting_character_field(index)
        return self._splitting_field(Keps)

    @cached_method
    def decomposition_field(self):
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

        OUTPUT:
        
        The composite field of :meth:`complete_definition_field` and
        :meth:`splitting_field`

        EXAMPLE::
        
            sage: K.<t> = QuadraticField(-3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.complete_definition_field()
            Number Field in tu0 with defining polynomial x^4 + 10*x^2 + 1
            sage: E.splitting_field()
            Number Field in tzeta0 with defining polynomial x^4 + 36
            sage: E.decomposition_field()
            Number Field in tu0tzeta0 with defining polynomial x^8 + 32*x^6 + 456*x^4 - 1408*x^2 + 10000

        """
        Ksplit = self.splitting_field()
        Kphi = self.complete_definition_field()
        K = self.definition_field()
        Kd = self.degree_field()
        to_Kphi = K.hom(Kphi) * Kd.embeddings(K)[0]
        to_Ksplit = Kd.embeddings(Ksplit)[0]
        _, _, from_Kphi = composite_field(to_Ksplit, to_Kphi, give_maps=True)
        return write_as_extension(from_Kphi)

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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.does_decompose()
            False
            sage: E.decomposable_twist().does_decompose()
            True

        """
        if not self._is_cached("_beta"):
            self.splitting_map(verbose=-1)
        c = self.c
        c_beta = self.c_splitting_map
        G = self.decomposition_field().galois_group()
        return all(c(s,t) == c_beta(s,t) for s in G for t in G)

    def _splitting_map_first_guess(self):
        r"""Give a naive guess of a splitting map.

        """
        eps = self.splitting_character(galois=True)
        d = self.degree_map
        Lbeta = self.splitting_image_field()
        @cached_function(key=lambda s: (str(s), s.parent().number_field()))
        def beta(sigma):
            return sqrt(Lbeta(d(sigma) * eps(sigma)))
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
            if a == 1:
                return [0]
            elif a == -1:
                return [1]
            else:
                raise ValueError("%s is not 1 or -1"%a)
        try:
            alpha = function_with_coboundary(G, (1, [-1], [2], convert), c_err)
            beta0 = self._beta
            @cached_function(key=lambda s: (str(s), s.parent().number_field()))
            def beta(sigma):
                return beta0(sigma) * alpha(sigma)
            self._beta = beta
            self.c_splitting_map.clear_cache() # Delete values of previous beta
            if not self.does_decompose():
                raise ValueError("Should be impossible to reach this code!");
        except ArithmeticError:
            if verbose >= 0:
                print ("Warning: The restriction of scalars of this Q-curve " +
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
            return Lbeta(from_Lbeta0(beta0(sigma)) * from_Lchi(chi(sigma)))
        return beta

    @cached_method(key=lambda self, i, v: i)
    def splitting_map(self, index=0, verbose=False):
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
            [ 1 -2  1  2 -2  1  2  1]
            [ 1  1 -1 -1  1 -1 -1  1]
            [ 1  2 -1  2  2 -1  2  1]
            [ 1 -2  1  2 -2  1  2  1]
            [ 1  1 -1 -1  1 -1 -1  1]
            [ 1  2 -1  2  2 -1  2  1]
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
        if hasattr(index, "__iter__"):
            return tuple(self.splitting_map(i, verbose=verbose) for i in index)
        if index in ZZ:
            if index == 0:
                return self._first_splitting_map(verbose=verbose)
            else:
                return self._indexed_splitting_map(index)
        if index == 'all':
            n = self.number_of_splitting_maps()
            return self.splitting_map(tuple(range(n)), verbose=verbose)
        if index == 'conjugacy':
            return tuple(self.splitting_map(index=ii[0], verbose=verbose)
                         for ii in self._conjugacy_classes())
        raise Exception("Invalid index %s"%index)

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
            [ 1 -2  1  2 -2  1  2  1]
            [ 1  1 -1 -1  1 -1 -1  1]
            [ 1  2 -1  2  2 -1  2  1]
            [ 1 -2  1  2 -2  1  2  1]
            [ 1  1 -1 -1  1 -1 -1  1]
            [ 1  2 -1  2  2 -1  2  1]
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
            self.splitting_map(verbose=verbose);
        return QQ(self._beta(sigma) * self._beta(tau) *
                  self._beta(sigma*tau)^(-1))

    def _Kphi_roots(self):
        r"""Give a basis for the field of complete definition.

        OUTPUT:
        
        A list of squarefree integers such that the field of complete
        definition is $\Q$ adjoint all roots of these
        integers. Furthermore this list has minimal length in this
        regard.

        """
        Kphi = self.complete_definition_field()
        products = [1]
        result = []
        for tmp in Kphi.subfields(degree=2):
            c = tmp[0].discriminant().squarefree_part()
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
        decomposition field.

        .. NOTE::

        For this $N$ to exist, the definition field of this Q-curve
        should be abelian.

        .. SEEALSO::

            :meth:`decomposition_field`

        OUTPUT:
        
        The smallest non-negative integer $N$ such that the
        decomposition field of this Q-curve as given by
        :meth:`decomposition_field` is a subfield of $\Q(\zeta_N)$.

        EXAMPLE::

            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: N = E.cyclotomic_order(); N
            40
            sage: L = E.decomposition_field(); L
            Number Field in tu0zeta0 with defining polynomial x^8 - 22*x^6 - 20*x^5 + 199*x^4 + 380*x^3 + 882*x^2 + 740*x + 1721
            sage: len(L.gen().minpoly().change_ring(CyclotomicField(N)).roots()) > 0
            True

        """
        if not self._is_cached('_N') or not self._is_cached('_ker'):
            K = self.decomposition_field()
            N = K.conductor(check_abelian=True)
            self._N = N
            L.<zeta> = CyclotomicField(N)
            a = K.gen().minpoly().change_ring(L).roots()[0][0]
            self._ker = [n for n in range(N) if gcd(n, N) == 1 and
                         a == sum(ai * zeta^(i*n)
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
                if chi(x) != 1:
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
        if hasattr(i, "__iter__"):
            return [self._twist_character_data(ii, j) for ii in i]
        if i in ZZ:
            if j == 1 and len(self._chi[i]) < 2:
                self._chi[i].append(dirichlet_to_galois(self._chi[i][0]))
            return self._chi[i][j]
        if i == 'all':
            return tuple(self._twist_character_data(ii, j)
                         for ii in range(self.number_of_splitting_maps()))
        if i == 'conjugacy':
            return tuple(self._twist_character_data(ii[0], j)
                         for ii in self._conjugacy_classes())
        raise Exception("Invalid index %s."%i)
    
    def twist_character(self, index=0, galois=False):
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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12*t - 12, 0, 180 - 108*t, 0], guessed_degrees=[2])
            sage: beta = E.splitting_map('all')
            sage: xi = E.twist_character('all', galois=True)
            sage: G = E.decomposition_field().galois_group()
            sage: all(beta[i](s) == xi[i](s) * beta[0](s) for s in G for i in range(len(beta)))
            True

        """
        if galois:
            return self._twist_character_data(index, 1)
        else:
            return self._twist_character_data(index, 0)

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
        while len(beta_del) > 0:
            beta0 = beta_del[0]
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
                            result[-1].append(beta_dict[beta])
                            break
            for j in result[-1]:
                beta_del.remove(beta_ls[j])
        return tuple(result)

    def _isogeny_data(self, iota):
        r"""Give the isogeny data of this curve over a given field.

        INPUT:

        - ``iota`` -- A field homomorphism from the definition field
          of this Q-curve to another number field.

        OUTPUT:

        A dictionary of which the keys are the elements of the galois
        group of K, and the value for each element sigma is a tuple of
        a rational function in $x$ and an algebraic number, that are
        respectively the $x$-coordinate map of the isogeny from the
        sigma conjugate of this curve to this curve and the induced
        scalar multiplication on the differentials.

        """
        G = K.galois_group()
        F = self.isogeny_x_map
        l = self.isogeny_scalar
        Kphi = self.complete_definition_field()
        _, to_L, from_Kphi = composite_field(iota, Kphi,
                                             give_maps=True)
        L, to_L = write_as_extension(to_L, give_map=True)
        from_Kphi = to_L * from_Kphi
        return {s : (F(s).change_ring(from_Kphi), from_Kphi(l(s)))
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
                if len(K_to_L) > 0:
                    K_to_L = K_to_L[0]
                else:
                    K_to_L = None
            if K_to_L != None:
                return Qcurve(result, isogenies=self._isogeny_data(K_to_L))
        return result

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

            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2]); E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 5
            sage: E2 = E.twist(t); E2
            Q-curve defined by y^2 = x^3 + (1/4*sqrt_a0^3-3*sqrt_a0^2-11*sqrt_a0+76)*x^2 + (15/8*sqrt_a0^3-45/2*sqrt_a0^2-165/2*sqrt_a0+660)*x over Number Field in sqrt_a0 with defining polynomial x^4 - 16*x^3 - 8*x^2 + 576*x - 1264
            sage: E2.complete_definition_field()
            Number Field in sqrt_a0 with defining polynomial x^4 - 16*x^3 - 8*x^2 + 576*x - 1264

        """
        K_gamma = gamma.parent()
        Kdef = self.definition_field()
        Kphi = self.complete_definition_field()
        K, iota, gamma_map = composite_field(Kphi, K_gamma,
                                             give_maps=True)
        gamma = gamma_map(gamma)
        E_map = iota * Kdef.hom(Kphi)
        E = twist_elliptic_curve(self.change_ring(E_map), gamma)
        ainvs = E.a_invariants()
        R.<x> = K[]
        G = K.galois_group()
        l = self.isogeny_scalar
        F = self.isogeny_x_map
        isogenies=dict()
        for s in G:
            Fs = F(s).change_ring(iota)
            ls = l(s).change_ring(iota)
            x, y = Fs.parent().gens()
            Ls.<agamma> = K.extension(x^2 - s(gamma)/gamma)
            Fs = gamma * Fs(x / s(gamma), y)
            ls = agamma * ls
            isogenies[s] = (Fs.change_ring(Ls), ls)
        # TODO: minimization of resulting field
        return Qcurve(ainvs, isogenies=isogenies)

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
            CI = product(Pgen[i]^k[i] for i in range(len(k)))
            if CI not in H:
                H.append(CI)
        S0 = []
        skip = copy(H)
        for CI in CG:
            if CI not in skip and CI^2 in H:
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

            sage: K.<t> = QuadraticField(7)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E.does_decompose()
            False
            sage: E.decomposable_twist().does_decompose()
            True

        """
        K = self.decomposition_field()
        if self.does_decompose():
            return K(1)
        S = self._decomposable_twist_set()
        G = K.galois_group()
        US = K.S_unit_group(proof=False, S=S)
        def c_err(sigma, tau):
            return US(self.c(sigma, tau) / self.c_splitting_map(sigma, tau))
        alpha = function_with_coboundary(G, US, c_err)
        return hilbert90(K, lambda s: alpha(s)^2)
    
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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2]); E
            Q-curve defined by y^2 = x^3 + 12*x^2 + (18*t+18)*x over Number Field in t with defining polynomial x^2 - 3
            sage: E.degree_map_image()
            [1, 2]
            sage: K1 = E.complete_definition_field(); K1
            Number Field in tu0 with defining polynomial x^4 - 2*x^2 + 25
            sage: K1(2).is_square()
            False
            sage: E2 = E.complete_definition_twist([2]); E2
            Q-curve defined by y^2 = x^3 + (-12*t)*x^2 + (54*t+54)*x over Number Field in t with defining polynomial x^2 - 3
            sage: K2 = E2.complete_definition_field(); K2
            Number Field in tsqrt_a0 with defining polynomial x^4 - 10*x^2 + 1
            sage: K2(2).is_square()
            True

        """
        # Calculate all elements generated by absolute roots mod squares
        # and how to obtain them
        roots_image = {1 : []}
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
        # TODO update this function
        Kbase = self.definition_field()
        Kold = self.complete_definition_field()
        Kroots = QQ
        for i, a in enumerate(roots):
            Kroots = field_with_root(Kroots, a, names=('sqrt_a'+str(i)))
        base_to_old = self._to_Kl
        Knew, base_to_new, roots_to_new = composite_field(Kbase, Kroots,
                                                          give_maps=True)
        Kbig, old_to_big, new_to_big = composite_field(Kold, Knew,
                                                       give_maps=True)
        base_to_big = new_to_big * base_to_new

        # The map we want as scalars for the new curve
        d = self.degree_map
        @cached_function(key=lambda s: (str(s), s.parent().number_field()))
        def mu(s):
            my_iter = (roots[i] for i in roots_image[d(s).squarefree_part()])
            return sqrt(Knew(product(my_iter)))

        # The correction map
        l = self.isogeny_scalar
        def alpha(s):
            return new_to_big(mu(s))^2 / old_to_big(l(s))^2

        # The twist parameter
        gamma = hilbert90(Kbig, alpha)

        # Check if we can twist by an element of Kbase
        gamma_ls = [x for x,e in gamma.minpoly().change_ring(Kbase).roots()
                    if base_to_big(x) == gamma]
        if len(gamma_ls) > 0:
            gamma = gamma_ls[0]
            isogenies = {s : (mu(s), d(s)) for s in Kbase.galois_group()}
            return Qcurve(twist_elliptic_curve(self, gamma),
                          isogenies=isogenies)

        # General case
        print ("Warning: Chosen twist is not defined over the same field " +
               "anymore.")
        E = twist_elliptic_curve(self.change_ring(base_to_big), gamma)
        ainvs = [[x for x,e in a.minpoly().change_ring(Knew).roots()
                  if new_to_big(x) == a] for a in E.a_invariants()]
        if product(len(a) for a in ainvs) == 0:
            raise ArithmeticError("The sought twist is not defined over " +
                                  "the given field.")
        ainvs = [a[0] for a in ainvs]
        isogenies = {s : (mu(s),d(s)) for s in Kbase.galois_group()}
        return Qcurve(ainvs, isogenies=isogenies)

    @cached_method
    def conductor_restriction_of_scalars(self):
        r"""Give the conductor of the restriction of scalars of this Q-curve.

        OUTPUT:

        The conductor of the restriction of scalars of this curve over
        the decomposition field.

        EXAMPLE::

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
            iota = K0.hom(K)
            # Proposition 1 of Milne, On the arithmetic of Abelian varieties
            return (self.change_ring(iota).conductor().absolute_norm() *
                    K.discriminant()^2)
        else:
            return self.conductor().absolute_norm() * K.discriminant()^2

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

        """
        if N is None:
            N = self.conductor_restriction_of_scalars()
        eps = [c^(-1) for c in self.splitting_character(index='conjugacy')]
        chi = [c^(-1) for c in self.twist_character(index='conjugacy')]
        M = lcm(character.modulus() for character in (eps + chi))
        L = QQ
        for character in eps + chi:
            L = composite_field(L, character.base_ring())
        D = DirichletGroup(M, base_ring=L)
        eps = [D(character) for character in eps]
        chi = [D(character) for character in chi]
        alpha = tuple(eps[i].conductor() for i in range(len(eps)))
        beta = tuple(tuple((chi[j] * chi[i]^(-1)).conductor()
                           for j in range(len(eps))) for i in range(len(eps)))
        gamma = tuple(tuple((chi[j] * chi[i]^(-1) * eps[i]).conductor()
                            for j in range(len(eps))) for i in range(len(eps)))
        d = [Kf.degree()
             for Kf in self.splitting_image_field(index='conjugacy')]
        primes = ZZ(N).prime_factors()
        level_dict = {p : self._newform_levels(prime=p, alpha=alpha, beta=beta,
                                               gamma=gamma, d=d, N=N)
                      for p in primes}
        return [tuple(product(primes[i]^level_dict[primes[i]][x[i]][j]
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
        x_max = max(max(max(beta[i][j] + 1, beta[i][j] + gamma[i][j])
                        for j in range(len(beta[i])))
                    for i in range(len(beta)))
        x_ls = []
        if sum(d) * x_max >= N: # Only small possibilities
            for x in mrange([x_max + 1]*len(alpha)): 
                candidate = (sum(d[i] * x[i] for i in range(len(d))) == N)
                for i in range(len(beta)):
                    if not candidate:
                        break
                    candidate = (alpha[i] <= x[i])
                    for j in range(len(beta[i])):
                        if not candidate:
                            break
                        if beta[i][j] == 0:
                            candidate = (x[j] == x[i])
                        elif (x[i] > max(beta[i][j] + 1,
                                         beta[i][j] + gamma[i][j]) or
                              (x[i] < max(beta[i][j] + 1,
                                          beta[i][j] + gamma[i][j]) and
                               gamma[i][j] >= 2) or
                              (x[i] == 1 and alpha[i] == 1 and
                               beta[i][j] == 1 and gamma[i][j] == 1)):
                            candidate = (x[j] == max(x[i], beta[i][j] + 1,
                                                     beta[i][j] + gamma[i][j]))
                        else:
                            candidate = (x[j] <= max(x[i], beta[i][j] + 1,
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

    def trace_of_frobenius(self, prime, power=1, splitting_map=0):
        r"""Compute the trace of a Frobenius element acting on this curve.

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
        not ramify in the decomposition field of this curve; this
        curve has good reduction at some prime $P$ above $p$ of the
        decomposition field; and the reduction of the dual of $\phi$
        to this good reduction is separable, where $\phi$ is the
        isogeny of the conjugate of this curve by a $p$-Frobenius
        element to this curve. In this case the corresponding
        $\lambda$-adic representation is unramified at
        $p$. Furthermore the trace of the Frobenius element at $p$ can
        be explicitly computed and is an algebraic number that only
        depends on $\beta$, the reduction at $P$ and the dual isogeny
        of $\phi$.

        This function computes the algebraic number that is the trace
        of (a power of) a Frobenius element. It will first check that
        the give prime number `prime` satisfies the mentioned
        conditions and raise a ValueError if this is not the case.

        INPUT:

        - ``prime`` -- A prime number such that `prime` does not
          ramify in the decomposition field of this curve; this curve
          has good reduction at some prime $P$ above `prime` in the
          decomposition field; and the reduction of the dual $\phi$ to
          this good reduction is separable, where $\phi$ is the
          isogeny from the conjugate of this curve by a Frobenius
          element of `prime` to this curve itself.

        - ``power`` -- A strictly positive integer (default: 1). If
          set to a value higher than 1 will compute the trace of the
          Frobenius element to the given power instead of the
          Frobenius element itself.

        - ``splitting_map`` -- A non-negative integer smaller than the
          number of splitting maps associated to this Q-curve
          (default: 0). This indicates the splitting map for which the
          trace of frobenius should be computed.

        OUTPUT:

        An algebraic number that is the trace of the Frobenius element
        at `prime` to the power `power` under any $\lambda$-adic
        Galois representation associated to the splitting map given by
        `splitting_map` for $\lambda$ not dividing `prime`.

        """
        if power > 1:
            pass # TODO
        else:
            K = self.decomposition_field()
            iota = self.definition_field().hom(K)
            E = self.change_ring(iota)
            if prime.divides(K.absolute_discriminant()):
                raise ValueError("The decomposition field is ramified " +
                                 "at " + str(prime))
            for P in K.primes_above(prime):
                if E.has_good_reduction(P):
                    break
            else:
                raise ValueError("This curve does not have good " +
                                 "reduction at any prime above " +
                                 str(prime))
            #TODO

    @cached_method
    def newform(self, algorithm='magma', verify=0):
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
        newforms of the lowest level is computed with the appropiate
        character. For these newforms the L-series of this Q-curve and
        the L-series of the newforms are compared until as many
        newforms are left as expected. Some Q-curves have multiple
        factors of the same level and character, so one is chosen.

        INPUT:
        
        - ``algorithm`` -- A string that determines which program
          should be used to do computations. Allowed values are:
          'magma' (default) or 'sage' to use MAGMA and Sage
          respectively.  This program is used for computing newforms
          and Euler factors of the L-series associated to the curve.

        - ``verify`` -- A non-negative integer determining what the
          biggest prime is for which the result should be verified
          using the Euler factors of the L-series.

        OUTPUT:

        A tuple consisting of

        - a newform such that this curve is a quotient of the abelian
          variety associated to that newform, i.e. a newform
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

            sage: K.<t> = QuadraticField(3)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E2 = E.decomposable_twist()
            sage: E2.newform() # long
            (q + (-a + 1)*q^3 + a*q^5 + 3*a*q^7 + (-2*a - 1)*q^9 + 4*q^11 + O(q^12),
             [Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1])

        Or another example with two factors::

            sage: K.<t> = QuadraticField(5)
            sage: E = Qcurve([0, 12, 0, 18*(1 + t), 0], guessed_degrees=[2])
            sage: E2 = E.decomposable_twist()
            sage: E2.newform() # long
            (q + (a + 1)*q^7 - 4*a*q^11 + O(q^12),
             [Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1,
              Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -zeta4])

        """
        if not self.does_decompose():
            raise ValueError("Can not compute newform if the restriction of " +
                             "scalars does not decompose.")

        # Getting this curve over the decomposition field as E
        Kl = self.complete_definition_field()
        to_Kl = self._to_Kl
        K = self.decomposition_field()
        iota = Kl.embeddings(K)[0]
        E = self.change_ring(iota * to_Kl)
        
        levels = self.newform_levels()
        twists_base = self.twist_character('conjugacy')
        # Find a common base for all twists
        M = lcm(chi.modulus() for chi in twists_base)
        twists_base = [chi.extend(M)^(-1) for chi in twists_base]
        # Characters of the newforms
        eps_ls = [(eps^(-1)).primitive_character()
                  for eps in self.splitting_character('conjugacy')]
        Lbeta_ls = self.splitting_image_field('conjugacy') # coefficient fields

        use_magma = (algorithm == 'magma')
        if not use_magma and algorithm != 'sage':
            raise ValueError("%s is not a valid algorithm to use."%algorithm)

        candidates = []
        max_level = lcm(lcm(tmp) for tmp in levels)
        # Keeps track of the lcm of all N considered
        # the primes in these will be excluded in checking
        # the Euler factors.
        done_cases = []
        for k in range(len(levels)):
            # Newform with smallest level:
            i_min, N = min(enumerate(levels[k]), key=(lambda x: x[1]))
            chi = twists_base[i_min]
            # Twists relative to i_min
            twists = [chi_j * chi^(-1) for chi_j in twists_base]
            eps = eps_ls[i_min]
            Lbeta = Lbeta_ls[i_min]
            expected_matches = sum(eps_ls[i] in eps.galois_orbit() and
                                   levels[k][i] == N
                                   for i in range(len(levels[k])))
            # Prevent getting the same newform twice:
            if (N, eps.galois_orbit(), Lbeta) in done_cases: 
                continue
            done_cases.append((N, eps.galois_orbit(), Lbeta))

            if use_magma:
                Dm = magma.DirichletGroup(eps.conductor(),
                                          magma(eps.base_ring()))
                for eps_m in Dm.Elements():
                    candidate = True
                    for i in Integers(eps.conductor()).unit_gens():
                        i = ZZ(i)
                        candidate = (eps(i) == eps_m(i))
                        if not candidate: # Not the right one
                            break
                    if candidate: # Found the right one
                        break
                if not candidate:
                    raise ValueError("No matching magma dirichlet character " +
                                     "for %s"%(candidate,))
                # The right level:
                eps_m = magma.DirichletGroup(N, magma(eps.base_ring()))(eps_m)
                cfs = magma.CuspForms(eps_m)
                nfs = magma.Newforms(cfs)
                nfs = [f[1] for f in nfs]
                E = magma(E)
            else:
                nfs = Newforms(eps.extend(N), names='a')

            # See if the coefficient field matches good enough:
            for f in nfs:
                if f.parent() == magma:
                    Kf = f.BaseField().sage()
                else:
                    Kf = f.base_ring()
                if (Kf.absolute_degree() == Lbeta.absolute_degree() and
                    Kf.absolute_discriminant() == Lbeta.absolute_discriminant()):
                    candidates.append((f, twists, N, eps))

        valid_options = []
        for k in range(len(levels)):
            option = {}
            for i in range(len(levels[k])):
                key = (levels[k][i], eps_ls[i].galois_orbit()[0])
                if key in option:
                    option[key] += 1
                else:
                    option[key] = 1
            valid_options.append(option)

        def is_submultiset(small, big):
            for key in small:
                if key not in big or small[key] > big[key]:
                    return False
            return True
        
        def valid_candidates(candidates):
            if len(candidates) > len(eps_ls):
                return False
            sub_option = {}
            for f, twists, N, eps in candidates:
                key = (N, eps.galois_orbit()[0])
                if key in sub_option:
                    sub_option[key] += 1
                else:
                    sub_option[key] = 1
            return any(is_submultiset(sub_option, option)
                       for option in valid_options)
                    
        p = 1
        while p.divides(max_level):
            p = next_prime(p)
        while p < verify or not valid_candidates(candidates):
            if use_magma:
                P_E = E.EulerFactor(p).sage()
            else:
                P_E = Euler_factor_elliptic_curve(E, p)
            remove = []
            for f, twists, N, eps in candidates:
                P_f = Euler_factor_modular_form(f, p, twists=twists)
                if P_f.parent() == magma:
                    P_f = P_f.sage()
                if P_f.list() != P_E.list():
                    remove.append((f, twists, N, eps))

            for tmp in remove:
                candidates.remove(tmp)

            p = next_prime(p)
            while p.divides(max_level):
                p = next_prime(p)

        return candidates[0][0], candidates[0][1]

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
        s += "over %s"%self.definition_field()
        return s
