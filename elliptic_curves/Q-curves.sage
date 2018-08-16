from sage.schemes.elliptic_curves.ell_number_field import EllipticCurve_number_field

def _lambda_of_isogeny(phi):
    """
    For an isogeny of the form
    ..MATH:

        (x,y) \mapsto (F(x), yF'(x) / \lambda)

    returns $\lambda$
    """
    Fx, Fy = phi.rational_maps()
    x,y = Fy.parent().gens()
    R = Fy.parent().base().base()
    u =  R((Fx.derivative(x) * y) / Fy)
    f = u.minpoly()
    return NumberField(f, names='l').gen()

def _lambda_of_isomorphism(E1, E2):
    """
    For two isomorphic curves $E_1$ and $E_2$,
    let $u$ be the unique scaling factor such that
    the isomorphism can be given by
    ..MATH:

       (x,y) \mapsto (u^2 x, u^3 y)

    then returns 1/u.
    In this case the equations are of the form
    ..MATH:

       E_1 : y^2 = x^3 + A x + B
       E_2 : y^2 = x^3 + u^4 A x + u^6 B
    """
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
        raise Exception("Could not compute u for some reason.")
    f = un.minpoly()
    f = f(f.variables()[0]^n)
    f = f.factor()[0][0]
    L.<u> = NumberField(f)
    return 1/u

class Qcurve(EllipticCurve_number_field):
    r"""
    A Q-curve over some number field
    """
    def _is_cached(self, var):
        return hasattr(self, var) and getattr(self, var) != None

    def __init__(self, curve, isogenies={}, guessed_degrees=[]):
        r"""
        Constructor of a Q-curve

        Will build all data associated to a Q-curve in the following way.
         - First of all will initialize the curve itself by using either
           a given curve or by data that can be turned into a curve using
           the function EllipticCurve. The resulting curve should be
           defined over some number field and will be redefined over the
           minimal galois extension of this field.
         - Next it will compute all the galois conjugates of itself with
           respect to the galois group of its base field.
         - Third it will fill in data about the isogenies from a galois
           conjugate to itself using the provided data given.
         - Next is will determine the remaining isogenies using the
           guessed degrees and combining previously found isogenies.
        
        If in this way it can not find isogenies from each galois conjugate
        to itself, the initialization will produce an error and the resulting
        object should not be used.

        INPUT:
        
         - ``curve`` -- An elliptic curve over some number field or any
                        input that would create such a curve when passed
                        to the constructor EllipticCurve. This curve will
                        be taken over the minimal galois extension of its
                        base field.
         - ``isogenies`` -- A dictionary (default: {}) with as keys elements
                            of the galois group of the base field of the
                            Q-curve and as values data of the corresponding
                            isogeny from the galois conjugate of this Q-curve
                            to itself. This data can be either an isogeny
                            as a Sage object or a tuple of an algebraic integer
                            (defined as an element of some number field) and
                            a strictly positive integer, which are respectively
                            the $\lambda$ such that the isogeny is
                            $z \mapsto \lambda z$ on the complex numbers and
                            the degree of the isogeny.
         - ``guessed_degrees`` -- A list (default: []) of strictly positive
                                  integers indicating possible degrees of
                                  isogenies from galois conjugates of this
                                  curve to itself.
        """
        self._init_curve(curve)
        self._init_isogenies()
        for sigma, phi in isogenies.iteritems():
            self._add_isogeny(sigma, phi)
        flag = self._fill_isogenies()
        for d in guessed_degrees:
            self._add_isogenies_of_degree(d)
            flag = self._fill_isogenies()
            if flag:
                break
        if not flag:
            raise ValueError("There is not sufficient isogeny information to make %s a Q-curve"%curve)

    def _init_curve(self, curve):
        if not isinstance(curve, EllipticCurve_number_field):
            curve = EllipticCurve(curve)
        if not isinstance(curve, EllipticCurve_number_field):
            raise ValueError("%s can not be a Q-curve"%curve)
        K = curve.base_ring()
        Kgal = K.galois_closure(names=K.variable_name() + 'g')
        G = Kgal.galois_group()
        iota = K.hom([a.minpoly().change_ring(Kgal).roots()[0][0] for a in K.gens()], Kgal)
        ainvs = [iota(a) for a in curve.a_invariants()]
        EllipticCurve_number_field.__init__(self, Kgal, ainvs)

    @cached_method
    def galois_conjugate(self, sigma):
        r"""
        Gives the galois conjugate of this curve.

        INPUT:

        - ``sigma`` -- A galois homomorphism of some number field

        OUTPUT:
        
        The galois conjugate of this curve by the galois homomorphism
        which extends to a common galois homomorphism over the algebraic
        closure of Q as sigma. This will be an elliptic curve and not
        returned as a Q-curve
        """
        sigma = galois_field_change(sigma, self.base_ring())
        return conjugate_curve(self, sigma)

    # Isogeny related stuff
    def _init_isogenies(self):
        self._l = dict()
        self._d = dict()
        self._Kl = self.base_ring()
        self._to_Kl = self._Kl.hom(self._Kl)
        e = self._Kl.galois_group().identity()
        self._l[e] = QQ(1)
        self._d[e] = 1

    def _add_isogeny(self, sigma, phi):
        if isinstance(phi, tuple):
            self._l[sigma], self._d[sigma] = phi
            self._update_isogeny_field()
        else:
            self._add_isogeny(sigma, (_lambda_of_isogeny(phi), phi.degree()))

    def _update_isogeny_field(self):
        G = list(self.base_ring().galois_group())
        for i in range(len(G)):
            if G[i] in self._l and self._l[G[i]] != None and self._l[G[i]].parent() != self._Kl:
                self._Kl, old_to_new, i_to_new = composite_field(self._Kl, self._l[G[i]].parent(), give_maps=True)
                self._to_Kl = old_to_new * self._to_Kl
                self._l[G[i]] = i_to_new(self._l[G[i]])
                for j in range(i):
                    if G[j] in self._l and self._l[G[j]] != None:
                        self._l[G[j]] = old_to_new(self._l[G[j]])
        if not self._Kl.is_galois():
            self._Kl, clos = self._Kl.galois_closure(names='al', map=True)
            self._to_Kl = clos * self._to_Kl
            for s in G:
                if s in self._l and self._l[s] != None:
                    self._l[s] = clos(self._l[s])

    def _fill_isogenies(self):
        G = self.base_ring().galois_group()
        Kl = self._Kl
        for s in G:
            for t in G:
                if (s*t not in self._l or self._l[s*t] == None) and \
                   (s in self._l and self._l[s] != None) and \
                   (t in self._l and self._l[t] != None):
                    tL = galois_field_extend(t, self._Kl, )
                    self._l[s*t] = tL(self._l[s]) * self._l[t]
                    self._d[s*t] = self._d[s] * self._d[t]
        for s in G:
            flag = s in self._l and self._l[s] != None
            if not flag:
                return flag
        return flag

    def _add_isogenies_of_degree(self, degree):
        G = self.base_ring().galois_group()
        fd = self.torsion_polynomial(degree)
        Kd, yotad = fd.splitting_field(names='a'+str(degree), map=True)
        Ed = self.change_ring(yotad)
        for g,e in fd.change_ring(yotad).factor():
            psi = Ed.isogeny(g)
            E_t = psi.codomain()
            j_t = E_t.j_invariant()
            for s in G:
                if yotad(self.galois_conjugate(s).j_invariant()) == j_t:
                    print "Degree %s isogeny found for"%degree, s
                    l1 = _lambda_of_isomorphism(self.galois_conjugate(s).change_ring(yotad),E_t)
                    l2 = _lambda_of_isogeny(psi.dual())
                    Kl, p1, p2 = composite_field(l1.parent(), l2.parent(), give_maps=True)
                    self._add_isogeny(s, (p1(l1) * p2(l2), psi.degree()))

    def isogeny_lambda(self, sigma):
        r"""
        Returns the $\lambda$ of the isogeny from the sigma conjugate 
        of this curve to this curve.

        The $\lambda$ of an isogeny is the constant such that
        the isogeny becomes $z \mapsto \lambda z$ on the complex
        numbers.

        INPUT:
        
        - ``sigma`` -- A galois homomorphism of a number field

        OUTPUT:

        The constant $\lambda$ such that the map $z \mapsto \lambda z$
        on the complex number defines an isogeny from a galois
        conjugate of this curve to itself. The galois conjugate is one
        obtained by conjugating with an extension of sigma.
        """
        if sigma not in self._l:
            self._l[sigma] = self._l[galois_field_change(sigma, self.base_ring())]
        return self._l[sigma]

    def complete_definition_field(self):
        r"""
        Gives the field over which the Q-curve is completely defined.

        OUTPUT:

        A number field over which both this elliptic curve and all
        isogenies from its galois conjugates to itself are defined.
        """
        return self._Kl

    def degree_map(self, sigma):
        r"""
        Gives the degree of an isogeny from the sigma galois conjugate
        of this curve to this curve.

        INPUT:

        - ``sigma`` -- A galois homomorphism of a number field

        OUTPUT:
        
        The degree of an isogeny from a galois conjugate of this curve
        to this curve itself. The galois conjugate is one obtained
        by conjugating with an extension of the given galois homomorphism
        sigma.
        """
        if sigma not in self._d:
            self._d[sigma] = self._d[galois_field_change(sigma, self.base_ring())]
        return self._d[sigma]

    @cached_method
    def degree_map_image(self):
        r"""
        Gives the image of the degree map in $\Q^*/(\Q^*)^2$

        OUTPUT:
        
        A list of squarefree integers such that each value of
        the degree map is a square times such an integer and
        all integers in this list differ a square from a value
        of the degree map.
        """
        result = []
        d = self.degree_map
        G = self.base_ring().galois_group()
        for s in G:
            val = d(s).squarefree_part()
            if val not in result:
                result.append(val)
        return result

    def degree_field(self):
        r"""
        Gives the fixed field of the degree map.

        OUTPUT:

        The biggest number field such that for each galois homomorphism
        that acts trivially on this field the degree map takes a value
        in $\Q^2$.
        """
        Kerd = []
        d = self.degree_map
        G = self.base_ring().galois_group()
        for s in G:
            if d(s).is_square():
                Kerd.append(s)
        return fixed_field(Kerd)

    def dual_basis(self, a1=None):
        r"""
        Gives a dual basis for the degree map.

        INPUT:

         - ``a1`` -- Optional parameter (default: None). If
                     set to a non-square integer which square
                     root is part of the degree field, will
                     ensure that this is the first entry of
                     the first list returned.

        OUTPUT:
        
        A tuple containing
         - A list of squarefree integers such that their
           square roots generate the degree field. This list
           is of minimal length with respect to such lists.
         - A list of non-negative integers of the same length
           as the first, such that the i-th entry differs
           precisely a square from the degree map at any
           galois homomorphism that only changes the sign
           of the square root of the i-th entry of the firs
           list and not of the others in that list.
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
    
    @cached_method
    def c(self, sigma, tau):
        r"""
        The value of the 2-cocycle $c: Gal(\bar{\Q}/\Q)^2 \to \Q^*$
        associated to a Q-curve by Quer.

        INPUT:

        - ``sigma`` -- A galois homomorphism over a number field
        - ``tau`` -- A galois homomorphism over a number field

        OUTPUT:

        The value
        ..MATH::

            \lambda_\sigma \cdot \sigma(\lambda_tau) \cdot lambda_{\sigma \tau}^{-1}

        where $\sigma$ and $\tau$ are extensions of sigma and tau to
        $\bar{\Q}$ respectively and where $\lambda_\sigma$ is the function
        :meth:`Qcurve.isogeny_lambda` at $\sigma$.
        """
        l = self.isogeny_lambda
        sigma = galois_field_change(sigma, self.complete_definition_field())
        tau = galois_field_change(tau, self.complete_definition_field())
        return QQ(l(sigma) * sigma(l(tau)) * l(sigma*tau)^(-1))

    def c_pm(self, sigma, tau):
        r"""
        The sign of the 2-cocycle $c: Gal(\bar{\Q}/\Q)^2 \to \Q^*$
        associated to a Q-curve by Quer.

        INPUT:

        - ``sigma`` -- A galois homomorphism over a number field
        - ``tau`` -- A galois homomorphism over a number field

        OUTPUT:

        The sign of
        ..MATH::

            \lambda_\sigma \cdot \sigma(\lambda_tau) \cdot lambda_{\sigma \tau}^{-1}

        where $\sigma$ and $\tau$ are extensions of sigma and tau to
        $\bar{\Q}$ respectively and where $\lambda_\sigma$ is the function
        :meth:`Qcurve.isogeny_lambda` at $\sigma$.
        """
        return sign(self.c(sigma,tau))

    def c_abs(self, sigma, tau):
        r"""
        The absolute value of the 2-cocycle $c: Gal(\bar{\Q}/\Q)^2 \to \Q^*$
        associated to a Q-curve by Quer.

        INPUT:

        - ``sigma`` -- A galois homomorphism over a number field
        - ``tau`` -- A galois homomorphism over a number field

        OUTPUT:

        The absolute value of
        ..MATH::

            \lambda_\sigma \cdot \sigma(\lambda_tau) \cdot lambda_{\sigma \tau}^{-1}

        where $\sigma$ and $\tau$ are extensions of sigma and tau to
        $\bar{\Q}$ respectively and where $\lambda_\sigma$ is the function
        :meth:`Qcurve.isogeny_lambda` at $\sigma$.
        """
        return abs(self.c(sigma,tau))

    @cached_method
    def xi_pm(self):
        r"""
        Returns the brauer group representation of the invariant $\xi_\pm$.

        OUTPUT:

        A list of tuples of integers, such that $\xi_\pm$ as an element
        of $Br_2(\Q)$ is the product of the quaternion algebras given
        by each of these tuples.
        """
        ai, di = self.dual_basis()
        return [(ai[i], di[i]) for i in range(len(ai))]

    @cached_method
    def _xi_pm_primes(self):
        r"""
        Gives the primes at which the $\xi_\pm$ might locally not be 1
        """
        result = lcm([lcm(h) for h in self.xi_pm()]).prime_factors()
        if 2 not in result:
            result.insert(0,2)
        return result

    @cached_method
    def xi_pm_local(self, p):
        r"""
        Gives a representative for $\xi_\pm$ in $Br_2(\Q_p)$.

        INPUT:

        - ``p`` -- A prime number.

        OUTPUT:
        
        +1 or -1 depending on whether the central simple algebra
        associated to $\xi_\pm$ over $\Q_p$ is split or non-split
        respectively.
        """
        if p not in self._xi_pm_primes():
            return 1
        else:
            return product([hilbert_symbol(ai,di,p) for (ai,di) in self.xi_pm()])

    def _first_splitting_character(self):
        N = 1
        eps_ls = [DirichletGroup(1)[0]]
        for p in self._xi_pm_primes():
            if self.xi_pm_local(p) == -1:
                if p == 2:
                    N *= 4
                    eps_ls.append(DirichletGroup(4).gen())
                else:
                    N *= p
                    eps_ls.append(DirichletGroup(p).gen())
        return product([eps_p.extend(N) for eps_p in eps_ls]).primitive_character()

    def _splitting_character_data(self, i, j):
        r"""
        Manages data related to splitting characters, i.e. for each splitting
        character we have a list containing the following entries
         0 - the splitting character as a dirichlet character
         1 - the fixed field of that splitting character
         2 - the splitting character as a galois character on its fixed field
        """
        if not self._is_cached('_eps') or 0 not in self._eps:
            self._eps = dict()
            self._eps[0] = [self._first_splitting_character()]
        if hasattr(i, "__iter__"):
            return tuple(self._splitting_character_data(ii, j) for ii in i)
        if i in ZZ:
            if i not in self._eps:
                eps0 = self._eps[0][0]
                chi = self.twist_character(i).primitive_character()
                N = lcm(eps0.conductor(),chi.conductor())
                self._eps[i] = [(eps0.extend(N) * chi.extend(N)^2).primitive_character()]
            if j >= 1 and len(self._eps[i]) < 2:
                self._eps[i].append(dirichlet_fixed_field(self._eps[i][0]))
            if j >= 2 and len(self._eps[i]) < 3:
                self._eps[i].append(dirichlet_to_galois(self._eps[i][0]))
                                                        
            return self._eps[i][j]
        if i == 'all':
            return tuple(self._splitting_character_data(ii, j) for ii in range(self.number_of_splitting_maps()))
        if i == 'conjugacy':
            return tuple(self._splitting_character_data(ii[0], j) for ii in self._conjugacy_determination())
        raise Exception("Invalid index %s."%i)
    
    def splitting_character(self, index=0, galois=False):
        r"""
        Gives a splitting character of this Q-curve.

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
                       character to return. Accepted values are
                       non-negative integers smaller than
                       the total amount of splitting maps
                       or one of the special values:
                        'all' : for a tuple of all splitting
                                characters
                        'conjugacy' : for a tuple of splitting
                                      characters, one for each
                                      conjugacy class of
                                      splitting maps.
                       Also accepts tuples of accepted values
                       including tuples themselves.
        - ``galois`` -- A boolean (default: False) indicating
                        whether the splitting characters
                        should be given as galois or dirichlet
                        characters.

        OUTPUT:

        Returns the splitting character of the given index, given
        as a galois character of galois is set to True or as a
        Dirichlet character otherwise. If the index was 'all'
        or 'conjugacy' will return a tuple of such characters
        corresponding to the corresponding tuple of indices. If
        the given index was a tuple will return a tuple of outputs
        on each entry of this tuple in the same order.
        """
        if galois:
            return self._splitting_character_data(index, 2)
        else:
            return self._splitting_character_data(index, 0)

    def splitting_character_field(self, index=0):
        r"""
        Gives the fixed field of a splitting character of this Q-curve.

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
                       character of which this should be the
                       fixed field. Accepted values are
                       non-negative integers smaller than
                       the total amount of splitting characters
                       or one of the special values:
                        'all' : corresponding to all splitting
                                characters
                        'conjugacy' : corresponding to splitting
                                      characters, one for each
                                      conjugacy class of
                                      splitting maps.
                       Also accepts tuples of accepted values
                       including tuples themselves.

        OUTPUT:

        Returns the fixed field of a splitting character of the
        given index. If the index was 'all' or 'conjugacy' will
        return a tuple of such fields corresponding to the corresponding
        tuple of indices. If the given index was a tuple will return a
        tuple of outputs on each entry of this tuple in the same order.
        """
        return self._splitting_character_data(index, 1)

    def _splitting_image_field(self, eps, Keps):
        if isinstance(eps, tuple):
            return tuple(self._splitting_image_field(eps[i], Keps[i]) for i in range(len(eps)))
        b = None
        if 2.divides(Keps.degree()):
            b = Keps.subfields(degree=2)[0][0].discriminant().squarefree_part()
        ai, di = self.dual_basis(a1=b)
        L = CyclotomicField(2*eps.order())
        for i in range(len(di)):
            Kdi = QuadraticField(di[i])
            if ai[i] == b:
                Lbig, L_to_Lbig, Kdi_to_Lbig = composite_field(L, Kdi, give_maps=True)
                alpha = L_to_Lbig(L.gen()) * Kdi_to_Lbig(Kdi.gen())
                L = Lbig.subfield(alpha)[0]
            else:
                L = composite_field(L, Kdi)
        return L

    @cached_method
    def splitting_image_field(self, index=0):
        r"""
        Gives the image field of a splitting map of this Q-curve.

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
                       map of the wanted splitting field. Accepted
                       values are non-negative integers smaller
                       than the total amount of splitting maps
                       or one of the special values:
                        'all' : for a tuple of all splitting
                                fields
                        'conjugacy' : for a tuple of splitting
                                      fields, one for each
                                      conjugacy class of
                                      splitting maps.
                       Also accepts tuples of accepted values
                       including tuples themselves.

        OUTPUT:

        Returns the splitting field of the splitting map of the
        given index. If the index was 'all' or 'conjugacy' will
        return a tuple of such characters corresponding to the
        corresponding tuple of indices. If the given index was a
        tuple will return a tuple of outputs on each entry of this 
        tuple in the same order.
        """
        eps = self.splitting_character(index)
        Keps = self.splitting_field(index)
        return self._splitting_image_field(eps, Keps)

    def _splitting_field(self, Keps):
        if isinstance(Keps, tuple):
            return tuple(self._splitting_field(Keps_i) for Keps_i in Keps)
        Kd = self.degree_field()
        return composite_field(Kd, Keps)

    @cached_method
    def splitting_field(self, index=0):
        r"""
        Gives a splitting field of this Q-curve.

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
                       map corresponding to this splitting
                       field. Accepted values are non-negative
                       integers smaller than the total amount
                       of splitting maps or one of the special
                       values:
                        'all' : for a tuple of all splitting
                                characters
                        'conjugacy' : for a tuple of splitting
                                      characters, one for each
                                      conjugacy class of
                                      splitting maps.
                       Also accepts tuples of accepted values
                       including tuples themselves.
        - ``galois`` -- A boolean (default: False) indicating
                        whether the splitting characters
                        should be given as galois or dirichlet
                        characters.

        OUTPUT:

        Returns the splitting field corresponding to the
        splitting map of the given index. If the index was 'all'
        or 'conjugacy' will return a tuple of such characters
        corresponding to the corresponding tuple of indices. If
        the given index was a tuple will return a tuple of outputs
        on each entry of this tuple in the same order.
        """
        Keps = self.splitting_character_field(index)
        return self._splitting_field(Keps)

    @cached_method
    def decomposition_field(self):
        r"""
        Gives the field over which the restriction of scalars of
        this Q-curve decomposes as a product of abelian varieties
        of GL_2-type.

        OUTPUT:
        
        The composite field of :meth:`Qcurve.complete_definition_field`
        and :meth:`Qcurve.splitting_field`
        """
        return composite_field(self.complete_definition_field(), self.splitting_field())

    def does_decompose(self):
        c = self.c
        c_beta = self.c_splitting_map
        G = self.decomposition_field().galois_group()
        for s in G:
            for t in G:
                if c(s,t) != c_beta(s,t):
                    return False
        return True

    def _splitting_map_first_guess(self):
        eps = self.splitting_character(galois=True)
        d = self.degree_map
        Lbeta = self.splitting_image_field()
        @cached_function
        def beta(sigma):
            return sqrt(Lbeta(d(sigma) * eps(sigma)))
        return beta

    def _first_splitting_map(self):
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
            @cached_function
            def beta(sigma):
                return beta0(sigma) * alpha(sigma)
            self._beta = beta
            for sigma in G:
                for tau in G:
                    if self.c(sigma, tau) != beta(sigma) * beta(tau) * beta(sigma*tau)^(-1):
                        raise ValueError("Should be impossible to reach this code!");
        except ArithmeticError:
            print "Warning: The restriction of scalars of this Q-curve over the "+\
                  "decomposition field does not decompose into abelian varieties"+\
                  " of GL_2-type. Pleas use the method decomposable_twist to "+\
                  "find a twist that does."
        return self._beta

    def _indexed_splitting_map(self, i):
        beta0 = self.splitting_map()
        Lbeta0 = self.splitting_image_field()
        chi = self.twist_character(i, galois=True)
        Lchi = self.twist_character(i).base_ring()
        L = composite_field(Lbeta0, Lchi)
        Lbeta = self.splitting_image_field(i)
        @cached_function
        def beta(sigma):
            return Lbeta(L(beta0(sigma)) * L(chi(sigma)))
        return beta

    @cached_method
    def splitting_map(self, index=0):
        r"""
        Gives a splitting map of this Q-curve.

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
                       map to return. Accepted values are
                       non-negative integers smaller than
                       the total amount of splitting maps
                       or one of the special values:
                        'all' : for a tuple of all splitting
                                maps
                        'conjugacy' : for a tuple of splitting
                                      maps, one for each
                                      conjugacy class of
                                      splitting maps.
                       Also accepts tuples of accepted values
                       including tuples themselves.

        OUTPUT:

        Returns the splitting map of the given index. If the
        index was 'all' or 'conjugacy' will return a tuple of
        such characters corresponding to the corresponding
        tuple of indices. If the given index was a tuple will
        return a tuple of outputs on each entry of this tuple
        in the same order.
        """
        if hasattr(index, "__iter__"):
            return tuple(self.splitting_map(i) for i in index)
        if index in ZZ:
            if index == 0:
                return self._first_splitting_map()
            else:
                return self._indexed_splitting_map(index)
        if index == 'all':
            return self.splitting_map(tuple(range(self.number_of_splitting_maps())))
        if index == 'conjugacy':
            return tuple(self.splitting_map(index=ii[0]) for ii in self._conjugacy_determination())
        raise Exception("Invalid index %s"%index)

    def c_splitting_map(self, sigma, tau):
        r"""
        Evaluates the coboundary of a splitting map of this Q-curve.

        Note that this is independent of the chosen splitting map.
       
        INPUT:
        
        - ``sigma`` -- A galois homomorphism of a number field.
        - ``tau`` -- A galois homomorphism of a number field.

        OUTPUT:
        
        The value
        ..MATH::

            \beta(\sigma) \cdot \beta(\tau) \cdot \beta(\sigma \tau)^{-1}

        for $\beta$ a splitting map of this Q-curve and $\sigma$ and
        $\tau$ galois extensions to $\bar{\Q}$ of sigma and tau respectively.
        """
        if not self._is_cached('_beta'):
            self.splitting_map();
        return QQ(self._beta(sigma) * self._beta(tau) * self._beta(sigma*tau)^(-1))

    def _Kl_roots(self):
        r"""
        Gives a basis for the field of complete definition.

        OUTPUT:
        
        Gives a list of non-square integers such that the field
        of complete definition $\Q$ adjoint all roots of these
        integers. Furthermore this list has minimal length in
        this regard.
        """
        Kl = self._Kl
        products = [1]
        result = []
        for tmp in Kl.subfields(degree=2):
            c = tmp[0].discriminant().squarefree_part()
            if c not in products:
                result.append(c)
                products.extend([(c*b).squarefree_part() for b in products])
        if len(products) < Kl.degree():
            raise ValueError("This Q-curve is not completely defined over a 2-...-2 extension. "+\
                             "This method only works when it is.")
        return result

    def cyclotomic_order(self):
        r"""
        The smallest $N$ such that $\Q(\zeta_N)$ contains the decomposition field.

        OUTPUT:
        
        The smallest non-negative integer $N$ such that the decomposition field
        of this Q-curve as given by :meth:`Qcurve.decomposition_field` is
        completely contained in $\Q(\zeta_N)$.
        """
        if not self._is_cached('_N') or not self._is_cached('_ker'):
            ai = self._Kl_roots()
            eps_ls = [self.splitting_character()]
            eps_ls.extend([character_for_root(a) for a in ai])
            N = lcm([eps.modulus() for eps in eps_ls])
            ker_ls = [set(eps.extend(N).kernel()) for eps in eps_ls]
            ker = ker_ls[0]
            for i in range(1,len(ker_ls)):
                ker = ker.intersection(ker_ls[i])
            self._N = N
            self._ker = list(ker)
        return self._N

    def _init_twist_characters(self):
        N = self.cyclotomic_order()
        ker = self._ker
        D = DirichletGroup(N)
        self._chi = []
        for chi in D:
            flag = True
            for x in ker:
                if chi(x) != 1:
                    flag = False
                    break
            if flag:
                self._chi.append([chi])

    def _twist_character_data(self, i, j):
        r"""
        Keeps track of the twist characters, storing for each splitting
        map the following data
         0) The dirichlet character which twists the default splitting
            map into this one.
         1) The same character as a galois character.
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
            return tuple(self._twist_character_data(ii, j) for ii in range(self.number_of_splitting_maps()))
        if i == 'conjugacy':
            return tuple(self._twist_character_data(ii[0], j) for ii in self._conjugacy_determination())
        raise Exception("Invalid index %s."%i)
    
    def twist_character(self, index=0, galois=False):
        r"""
        Gives the twist needed to obtain a certain splitting map.

        INPUT:

        - ``index`` -- The index (default: 0) of the splitting
                       map to which this twist should correspond.
                       Accepted values are non-negative integers 
                       smaller than the total amount of splitting
                       maps or one of the special values:
                        'all' : for a tuple of all splitting
                                characters
                        'conjugacy' : for a tuple of splitting
                                      characters, one for each
                                      conjugacy class of
                                      splitting maps.
                       Also accepts tuples of accepted values
                       including tuples themselves.
        - ``galois`` -- A boolean (default: False) indicating
                        whether the twist character should be
                        given as galois or dirichlet characters.

        OUTPUT:

        Returns the character such that twisting the default
        splitting map by this character gives the splitting
        map of the given index. This twist character isgiven
        as a galois character if galois is set to True or as a
        Dirichlet character otherwise. If the index was 'all'
        or 'conjugacy' will return a tuple of such characters
        corresponding to the corresponding tuple of indices. If
        the given index was a tuple will return a tuple of outputs
        on each entry of this tuple in the same order.
        """
        if galois:
            return self._twist_character_data(index, 1)
        else:
            return self._twist_character_data(index, 0)

    @cached_method
    def number_of_splitting_maps(self, count_conjugates=True):
        r"""
        Gives the number of splitting maps for this Q-curve.

        INPUT:

        - ``count_conjugates`` -- A boolean (default: True)
                                  indicating whether conjugate
                                  splitting maps should be
                                  counted seperately.

        OUTPUT:
        
        Gives the number of distinct splitting maps of this
        Q-curve defined over its decomposition field. If the
        flag count_conjugates is set to False, will return
        the number of conjugacy classes of such splitting
        maps.
        """
        if count_conjugates:
            if not self._is_cached('_chi'):
                self._init_twist_characters()
            return len(self._chi)
        else:
            return len(self._conjugacy_determination)

    @cached_method
    def _conjugacy_determination(self):
        r"""
        Gives a tuple of indices that contains for each class
        of splitting maps the index of exactly one element
        thereof.
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
                            if tau(L0_to_M(beta0(sigma))) != L_to_M(beta(sigma)):
                                flag = False
                                break
                        if flag:
                            result[-1].append(beta_dict[beta])
                            break
            for j in result[-1]:
                beta_del.remove(beta_ls[j])
        return tuple(result)

    def twist(self, gamma):
        r"""
        Gives the twist of this Q-curve by a given element gamma.

        INPUT:

        - ``gamma`` -- An element of a number field.

        OUTPUT:
        
        A Q-curve defined over the composite field of the field over which this
        Q-curve is completely defined and the parent of gamma, that is the
        twist of this Q-curve by gamma, i.e. if this Q-curve was given by
        
        ..MATH::

        E : y^2 = x^3 + a_2 x^2 + a_4 x + a_6

        the twisted Q-curve is given by

        ..MATH::
        
        E : y^2 = x^3 + \gamma a_2 x^2 + \gamma^2 a_4 x + \gamma^3 a_6
        
        """
        K_E = self.complete_definition_field()
        K_gamma = gamma.parent()
        K, iota, gamma_map = composite_field(K_E, K_gamma, give_maps=True)
        gamma = gamma_map(gamma)
        E_map = iota * self._to_Kl
        E = twist_elliptic_curve(self.change_ring(E_map), gamma)
        l = self.isogeny_lambda
        d = self.degree_map
        G = K.galois_group()
        isogenies = dict()
        for s in G:
            L, K_to_L, alpha = field_with_root(K, s(gamma)/gamma, give_embedding=True)
            isogenies[s] = (K_to_L(iota(l(s))) * alpha, d(s))
        return Qcurve(E, isogenies=isogenies)
    
    def decomposable_twist(self):
        r"""
        Gives another Q-curve which restriction of scalars over the decomposition
        field decomposes as a product of abelian varieties of GL_2-type.

        OUTPUT:
        
        A Qcurve which is a twist of this curve and has the same decomposition
        field. When taking the restriction of scalars of this curve over
        the decomposition field the resulting abelian variety is isognenous
        to a product of Q-simple, non-Q-isogenous abelian varieties of
        GL_2-type.
        """
        K = self.decomposition_field()
        CG = K.class_group()
        Pgen = [CG(product(K.primes_above(p))) for p in K.discriminant().prime_factors()]
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
        S = [P for I in S0 for P in I.prime_factors()]
        G = K.galois_group()
        US = K.S_unit_group(S=S)
        def c_err(sigma, tau):
            return US(self.c(sigma, tau) / self.c_splitting_map(sigma, tau))
        alpha = function_with_coboundary(G, US, c_err)
        gamma = 0
        while gamma == 0:
            x = K.random_element()
            gamma = sum(s(x) / K(alpha(s))^2 for s in G)

        # Minimizing gamma in some sense
        I = product(P^(floor(e/2)) for P, e in K.ideal(gamma).factor())
        J = (CG(I)^(-1)).ideal()
        gamma2 = ((J*I).gens_reduced()[0])^2  # Most of the square part of gamma
        gamma = gamma / gamma2
        gamma *= gamma.denominator() # Make kind of integral
        # Find biggest u s.t. u^k divides the n-k'th coefficient of the minimal polynomial of gamma
        cn = gamma.minpoly().list()
        cn.reverse() # Decreasing order
        n = lcm(i for i in range(1,len(cn)) if cn[i] != 0)
        un = gcd(cn[i]^(n/i) for i in range(1,len(cn)) if cn[i] != 0)
        u = product(x^floor(e/n) for x,e in un.factor())
        gamma = gamma / u # Minimal polynomial has smallest possible integers!

        # Updating alpha to fit the new gamma
        def alpha(s):
            s = galois_field_change(s, gamma.parent());
            return sqrt(s(gamma) / gamma)
        # Check to make sure everything is still good (computationally intensive, maybe remove?)
        def c_check(s, t):
            return QQ(alpha(s) * s(alpha(t)) / alpha(s*t) / c_err(s,t))
        def convert(a):
            if a == 1:
                return [0]
            elif a == -1:
                return [1]
            else:
                raise ValueError("%s is not 1 or -1"%a)
        try:
            alpha_check = function_with_coboundary(G, (1, [-1], [2], convert), c_check)
        except ArithmeticError:
            raise ValueError("Something went terribly wrong!")

        return self.twist(gamma)

    def complete_definition_twist(self, roots):
        r"""
        Gives a twist of this curve completely defined over a given field.

        INPUT:
        
        - ``roots`` -- A list of rational numbers satisfying the
          following property: There exists a set of generators of the
          image of the degree map in $\Q^*/(\Q^*)^2$ such that each
          element in the list is plus or minus an element in this set.
        
        OUTPUT:

        A Qcurve that is a twist of this curve and satisfies
         - It is defined over the same base field $K$
         - It is completely defined over $K$ adjoint all roots
           of the rationals given in roots.
        """
        # Calculate all elements generated by absolute roots mod squares and how to obtain them
        roots_image = {1 : []}
        for i in range(len(roots)):
            if abs(roots[i]).squarefree_part() not in roots_image:
                for b in list(roots_image):
                    roots_image[abs(roots[i]*b).squarefree_part()] = roots_image[b] + [i]

        # Check if roots is valid:
        d_image = self.degree_map_image()
        flag = (len(roots_image) != len(d_image)) # At least enough elements
        for a in roots: # Check each root plus or minus one associated to d
            if flag:
                break
            flag = abs(a) not in d_image
        if flag:
            raise ValueError("The set %s does not give a valid set of roots"%roots)

        # Let's compute the fields and corresponding embeddings
        Kbase = self.base_ring()
        Kold = self._Kl
        Kroots = QQ
        for a in roots:
            Kroots = field_with_root(Kroots, a)
        base_to_old = self._to_Kl
        Knew, base_to_new, roots_to_new = composite_field(Kbase, Kroots, give_maps=True)
        Kbig, old_to_big, new_to_big = composite_field(Kold, Knew, give_maps=True)
        base_to_big = new_to_big * base_to_new

        # The map we want as lambda for the new curve
        d = self.degree_map
        @cached_function
        def mu(s):
            return sqrt(Knew(product(roots[i] for i in roots_image[d(s).squarefree_part()])))

        # The correction map
        l = self.isogeny_lambda
        @cached_function
        def alpha(s):
            return new_to_big(mu(s))^2 / old_to_big(l(s))^2

        # The twist parameter
        gamma = 0
        while gamma == 0:
            x = Kbig.random_element()
            gamma = sum(s(x)/alpha(s) for s in Kbig.galois_group())
        gamma *= gamma.denominator() # Make element of O_K
        # Find biggest u s.t. u^k divides the n-k'th coefficient of the minimal polynomial of gamma
        cn = gamma.minpoly().list()
        cn.reverse() # Decreasing order
        n = lcm(i for i in range(1,len(cn)) if cn[i] != 0)
        un = gcd(cn[i]^(n/i) for i in range(1,len(cn)) if cn[i] != 0)
        u = product(x^floor(e/n) for x,e in un.factor())
        gamma = gamma / u # Minimal polynomial has smallest possible integers!

        # Check if we can twist by an element of Kbase
        gamma_ls = [x for x,e in gamma.minpoly().change_ring(Kbase).roots() if base_to_big(x) == gamma]
        if len(gamma_ls) > 0:
            gamma = gamma_ls[0]
            isogenies = {s : (mu(s), d(s)) for s in Kbase.galois_group()}
            return Qcurve(twist_elliptic_curve(self, gamma), isogenies=isogenies)

        # General case
        print "Warning: Chosen twist is not defined over the same field anymore."
        E = twist_elliptic_curve(self.change_ring(base_to_big), gamma)
        ainvs = [[x for x,e in a.minpoly().change_ring(Knew).roots() if new_to_big(x) == a] for a in E.a_invariants()]
        if product(len(a) for a in ainvs) == 0:
            raise ArithmeticError("The sought twist is not defined over the given field.")
        ainvs = [a[0] for a in ainvs]
        isogenies = {s : (mu(s),d(s)) for s in Kbase.galois_group()}
        return Qcurve(ainvs, isogenies=isogenies)

    @cached_method
    def conductor_restriction_of_scalars(self):
        r"""
        Gives the conductor of the restriction of scalars.

        OUTPUT:

        The conductor of the restriction of scalars of this curve
        over the decomposition field.
        """
        K0 = self.base_ring()
        K = self.decomposition_field()
        iota = K0.hom([a.minpoly().change_ring(K).roots()[0][0] for a in K0.gens()], K)
        # Proposition 1 of Milne, On the arithmetic of Abelian varieties
        return self.change_ring(iota).conductor().absolute_norm() * K.discriminant()^2
    
    @cached_method
    def level_table(self, prime=None, what='max'):
        eps_ls = [eps^(-1) for eps in self.splitting_character(index='conjugacy')]
        chi_ls = [chi^(-1) for chi in self.twist_character(index='conjugacy')]
        result = [0] * len(eps_ls)
        for i in range(len(eps_ls)):
            result[i] = [0] * len(chi_ls)
            eps = eps_ls[i]
            chi0 = chi_ls[i]
            for j in range(len(chi_ls)):
                N = lcm([eps.modulus(), chi_ls[j].modulus(), chi0.modulus()])
                chi = chi_ls[j].extend(N) * chi0.extend(N)^(-1)
                eps_chi = eps.extend(N) * chi
                if prime is None:
                    beta = chi.conductor()
                    gamma = eps_chi.conductor()
                    maxlevel = lcm(beta*product(beta.prime_factors()), beta*gamma)
                else:
                    beta = chi.conductor().ord(prime)
                    gamma = eps_chi.conductor().ord(prime)
                    maxlevel = max(beta + 1, beta + gamma)
                if what == 'max':
                    result[i][j] = maxlevel
                elif what == 'beta':
                    result[i][j] = beta
                elif what == 'gamma':
                    result[i][j] = gamma
        return result

    def _newform_levels(self, prime=None, alpha=None, beta=None, gamma=None, d=None, N=None):
        r"""
        Gives the possible levels of newforms associated to this Q-curve.

        INPUT:
        
        - ``prime`` -- A prime number or None (default: None) indicating the exponent
          of which prime the level should be computed or None for the level itself.
        - ``alpha`` -- A tuple of non-negative integers (default: None) containing the
          conductors of the characters associated to the newforms. Will be computed
          from the splitting characters up to conjugacy if set to None.
        - ``beta`` -- A tuple of tuples of non-negative integers (default: None)
          containing at index i, j the conductor of the twist that turns the i-th
          newform into the j-th newform. Will be computed from the twist characters
          up to conjugacy if set to None.
        - ``gamma`` -- A tuple of tuples of non-negative integers (default: None)
          containing at index i, j the conductor of the product of the twist that
          turns the i-th newform into the j-th newform and the character of the
          i-th newform. Will be computed from the twist characters and splitting
          characters up to conjugacy if set to None
        - ``d`` -- A tuple of non-negative integers (default: None) containing the
          respective degrees of the fields in which the newforms have their
          coefficients.
        - ``N`` -- A non-negative integer (default: None) giving the conductor of
          the restriction of scalars of this Q-curve over the decomposition field.
          Will be computed using the corresponding method if set to None.
        
        OUTPUT:

        A list of tuples, each tuple representing one of the options for the levels
        of the newforms associated to this Q-curve. If a prime was given, these
        tuples will contain the respective exponent of the given prime for each
        newform. If no prime was given, they will contain the respective level of
        each newform.
        """
        # Calculate missing stuff:
        if alpha is None or gamma is None:
            eps = [character^(-1) for character in self.splitting_character(index='conjugacy')]
        if beta is None or gamma is None:
            chi = [character^(-1) for character in self.twist_character(index='conjugacy')]
        if alpha is None or beta is None or gamma is None:
            M = lcm(character.modulus() for character in (eps + chi))
        if alpha is None or gamma is None:
            eps = [character.extend(M) for character in eps]
        if beta is None or gamma is None:
            chi = [character.extend(M) for character in chi]
        if alpha is None:
            alpha = tuple(eps[i].conductor() for i in range(len(eps)))
        if beta is None:
            beta = tuple(tuple((chi[j] * chi[i]^(-1)).conductor() for j in range(len(eps)))
                                                                  for i in range(len(eps)))
        if gamma is None:
            gamma = tuple(tuple((chi[j] * chi[i]^(-1) * eps[i]).conductor() for j in range(len(eps)))
                                                                            for i in range(len(eps)))
        if d is None:
            d = [Kf.degree() for Kf in self.splitting_image_field(index='conjugacy')]
        if N is None:
            N = self.conductor_restriction_of_scalars()

        if prime is None: # Level case
            level_dict = {}
            primes = ZZ(N).prime_factors()
            for p in primes:
                level_dict[p] = self._newform_levels(prime=p, alpha=alpha, beta=beta,
                                                     gamma=gamma, d=d, N=N)
            return [tuple(product(primes[i]^level_dict[primes[i]][x[i]][j] # Specific factor
                                  for i in range(len(primes))) # All primes
                                  for j in range(len(d))) # All entries
                                  for x in mrange([len(level_dict[p]) for p in primes])] # All options
        
        else: # prime case
            alpha = tuple(alpha[i].ord(prime) for i in range(len(alpha)))
            beta = tuple(tuple(beta[i][j].ord(prime) for j in range(len(beta[i]))) for i in range(len(beta)))
            gamma = tuple(tuple(gamma[i][j].ord(prime) for j in range(len(gamma[i]))) for i in range(len(gamma)))
            N = N.ord(prime)
            # Small cases
            x_max = max(max(max(beta[i][j] + 1, beta[i][j] + gamma[i][j]) for j in range(len(beta[i])))
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
                            elif x[i] > max(beta[i][j] + 1, beta[i][j] + gamma[i][j]) \
                                 or (x[i] < max(beta[i][j] + 1, beta[i][j] + gamma[i][j]) and \
                                     gamma[i][j] >= 2) \
                                     or (x[i] == 1 and \
                                         alpha[i] == 1 and \
                                         beta[i][j] == 1 and \
                                         gamma[i][j] == 1):
                                candidate = (x[j] == max(x[i], beta[i][j] + 1, beta[i][j] + gamma[i][j]))
                            else:
                                candidate = (x[j] <= max(x[i], beta[i][j] + 1, beta[i][j] + gamma[i][j]))
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
