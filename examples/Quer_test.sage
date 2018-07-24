### Index looper
def index_loop(n):
    j = [0]*len(n)
    loop = True
    while loop:
        yield j
        i = 0
        j[i] += 1
        while j[i] >= n[i]:
            j[i] = 0
            i += 1
            if i >= len(n):
                loop = False
                break
            j[i] += 1
    raise StopIteration()

### The main class
class Quer_invariants(SageObject):
    def __init__(self, E):
        K0 = E.base_ring()
        Kgal = K0.galois_closure(names=K0.variable_name() + 'g')
        G0 = Kgal.galois_group()
        iota = K0.hom([a.minpoly().change_ring(Kgal).roots()[0][0] for a in K0.gens()], Kgal)
        a_ls = [iota(a) for a in E.a_invariants()]
        KH = fixed_field([sigma for sigma in G0 if product(sigma(a) == a for a in a_ls)])
        a_ls = [KH(a) for a in a_ls]
        self._E0 = EllipticCurve(a_ls)

    def _is_cached(self, var):
        return hasattr(self, var) and getattr(self, var) != None

    def elliptic_curve(self):
        return self._E0

    def galois_conjugates(self):
        if not self._is_cached('_E'):
            self._E = dict()
            G = self._E0.base_ring().galois_group()
            for sigma in G:
                self._E[sigma] = conjugate_curve(self._E0, sigma)
        return self._E

    def _init_isogenies(self):
        flag = not self._is_cached('_l') or not self._is_cached('_d') or \
               not self._is_cached('_Kl') or not self._is_cached('_to_Kl')
        if flag:
            self._l = dict()
            self._d = dict()
            self._Kl = self._E0.base_ring()
            self._to_Kl = self._Kl.hom(self._Kl)
            e = self._Kl.galois_group().identity()
            self._l[e] = QQ(1)
            self._d[e] = 1
        return flag
    
    def lambda_of_isogeny(self, phi):
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

    def lambda_of_isomorphism(self, E1, E2):
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

    def _update_isogeny_field(self):
        if not self._is_cached('_Kl'):
            self._Kl = self._E0.base_ring()
            self._to_Kl = self._Kl.hom(self._Kl)
        G = list(self._E0.base_ring().galois_group())
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

    def add_isogeny(self, sigma, degree, l):
        """Manually add an isogeny from sigma(E) to E of a given degree with corresponding lambda l"""
        self._init_isogenies()
        K0 = self._E0.base_ring()
        G0 = K0.galois_group()
        if sigma not in G0:
            if len(sigma.parent().number_field().gen().minpoly().change_ring(K0).roots()) > 0: # Is a subfield
                for tau in G0:
                    if galois_field_change(tau, sigma.parent().number_field()) == sigma:
                        self.add_isogeny(tau, degree, l)
            else:
                self.add_isogeny(galois_field_change(sigma, K0))
        else:
            self._l[sigma] = l
            self._d[sigma] = degree
            self._update_isogeny_field()

    def add_isogenies(self, degree):
        self._init_isogenies()
        G = self._E0.base_ring().galois_group()
        fd = self._E0.torsion_polynomial(degree)
        Kd, yotad = fd.splitting_field(names='a_d', map=True)
        Ed = self._E0.change_ring(yotad)
        E = self.galois_conjugates()
        for g,e in fd.change_ring(yotad).factor():
            psi = Ed.isogeny(g)
            E_t = psi.codomain()
            j_t = E_t.j_invariant()
            for s in G:
                if yotad(E[s].j_invariant()) == j_t:
                    print "Case found for", s
                    l1 = self.lambda_of_isomorphism(E[s].change_ring(yotad),E_t)
                    l2 = self.lambda_of_isogeny(psi.dual())
                    Kl, p1, p2 = composite_field(l1.parent(), l2.parent(), give_maps=True)
                    self.add_isogeny(s, psi.degree(), p1(l1) * p2(l2))
        
    def fill_isogenies(self):
        self._init_isogenies()
        K = self._E0.base_ring()
        Kl = self._Kl
        G = K.galois_group()
        for s in G:
            for t in G:
                if (s*t not in self._l or self._l[s*t] == None) and \
                   (s in self._l and self._l[s] != None) and \
                   (t in self._l and self._l[t] != None):
                    tL = galois_field_extend(t, Kl)
                    self._l[s*t] = tL(self._l[s]) * self._l[t]
                    self._d[s*t] = self._d[s] * self._d[t]

    @cached_method
    def isogeny_lambda(self, sigma):
        if not self._is_cached('_l') or sigma not in self._l or self._l[sigma] == None:
            self.fill_isogenies()
            K = self._E0.base_ring()
            sigma = galois_field_change(sigma, K)
            if sigma not in self._l or self._l[sigma] == None:
                raise Exception("Please add some isogenies by using 'add_isogenies' with a guessed degree")
        return self._l[sigma]

    @cached_method
    def degree_map(self, sigma):
        if not self._is_cached('_d') or sigma not in self._d or self._d[sigma] == None:
            self.fill_isogenies()
            K = self._E0.base_ring()
            sigma = galois_field_change(sigma, K)
            if sigma not in self._d or self._d[sigma] == None:
                raise Exception("Please add some isogenies by using 'add_isogenies' with a guessed degree")
        return self._d[sigma]

    @cached_method
    def c(self, sigma, tau):
        if not self._is_cached('_c'):
            self._c = dict()
        self.fill_isogenies()
        s = galois_field_change(sigma, self._Kl)
        t = galois_field_change(tau, self._Kl)
        stpair = (s, t)
        if stpair not in self._c or self._c[stpair] is None:
            l = self.isogeny_lambda
            self._c[stpair] = QQ(l(s) * s(l(t)) * l(s*t)^(-1))
        return self._c[stpair]

    @cached_method
    def degree_field(self):
        Kerd = []
        d = self.degree_map
        K = self._E0.base_ring()
        G = K.galois_group()
        for s in G:
            if d(s).is_square():
                Kerd.append(s)
        Kd = fixed_field(Kerd)
        # Kerd = G.subgroup(Kerd)
        # Kd = Kerd.fixed_field()
        if isinstance(Kd, tuple):
            return Kd[0]
        else:
            return Kd
        
    def dual_basis(self, a1=None):
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
    def brauer(self):
        ai, di = self.dual_basis()
        return [(ai[i], di[i]) for i in range(len(ai))]

    @cached_method
    def xi_pm_relevant_primes(self):
        result = lcm([lcm(h) for h in self.brauer()]).prime_factors()
        if 2 not in result:
            result.insert(0,2)
        return result

    @cached_method
    def xi_pm(self, p):
        if p not in self.xi_pm_relevant_primes():
            return 1
        else:
            return product([hilbert_symbol(ai,di,p) for (ai,di) in self.brauer()])

    def _first_splitting_character(self):
        N = 1
        eps_ls = []
        for p in self.xi_pm_relevant_primes():
            if self.xi_pm(p) == -1:
                if p == 2:
                    N *= 4
                    eps_ls.append(DirichletGroup(4).gen())
                else:
                    N *= p
                    eps_ls.append(DirichletGroup(p).gen())
        return product([eps_p.extend(N) for eps_p in eps_ls]).primitive_character()

    def _splitting_character_data(self, i, j):
        if not self._is_cached('_eps') or 0 not in self._eps:
            self._eps = dict()
            self._eps[0] = [self._first_splitting_character()]
        if hasattr(i, "__iter__"):
            return [self._splitting_character_data(ii, j) for ii in i]
        if i in ZZ:
            if i not in self._eps:
                eps0 = self._eps[0][0]
                chi = self.twist_character(0).primitive_character()
                N = lcm(eps0.conductor(),chi.conductor())
                self._eps[i] = [eps0.extend(N) * chi.extend(N)^2]
            if j >= 1 and len(self._eps[i]) < 2:
                self._eps[i].append(dirichlet_fixed_field(self._eps[i][0]))
            if j >= 2 and len(self._eps[i]) < 3:
                self._eps[i].append(dirichlet_to_galois(self._eps[i][0],
                                                        self._eps[i][1]))
            return self._eps[i][j]
        if i == 'all':
            return [self._splitting_character_data(ii, j) for ii in range(self.number_of_splitting_maps())]
        if i == 'conjugacy':
            return [self._splitting_character_data(ii[0], j) for ii in self._conjugacy_determination()]
        raise Exception("Invalid index %s."%i)
    
    def splitting_character(self, index=0, galois=False):
        if galois:
            return self._splitting_character_data(index, 2)
        else:
            return self._splitting_character_data(index, 0)

    def splitting_character_field(self, index=0):
        return self._splitting_character_data(index, 1)

    def _splitting_image_field(self, eps, Keps):
        if isinstance(eps, list):
            return [self._splitting_image_field(eps[i], Keps[i]) for i in range(len(eps))]
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
        eps = self.splitting_character(index)
        Keps = self.splitting_field(index)
        return self._splitting_image_field(eps, Keps)

    def _splitting_field(self, Keps):
        if isinstance(Keps, list):
            return [self._splitting_field(Keps_i) for Keps_i in Keps]
        Kd = self.degree_field()
        return composite_field(Kd, Keps)
    
    @cached_method
    def splitting_field(self, index=0):
        Keps = self.splitting_character_field(index)
        return self._splitting_field(Keps)

    def _splitting_map_first_guess(self):
        eps = self.splitting_character(galois=True)
        d = self.degree_map
        Lbeta = self.splitting_image_field()
        @cached_function
        def beta(sigma):
            return sqrt(Lbeta(d(sigma) * eps(sigma)))
        return beta

    def _map_coboundary(self, f):
        @cached_function
        def d_f(sigma, tau):
            return f(sigma) * f(tau) * f(sigma*tau)^(-1)
        return d_f
    
    def _first_splitting_map(self):
        beta_guess = self._splitting_map_first_guess()
        c_beta = self._map_coboundary(beta_guess)
        c = self.c
        K = self.complete_field()
        @cached_function
        def c_cor(sigma, tau):
            return QQ(c(sigma, tau) / c_beta(sigma, tau))

        # Computations on the galois group
        G = K.galois_group()
        sqr_gens = []
        sqrs = []
        for sigma in G:
            if sigma^2 not in sqrs:
                sqrs.append(sigma^2)
                sqr_gens.append(sigma)
        mod_sqr_gens = []
        elim = copy(sqrs)
        for sigma in G:
            if sigma not in elim:
                mod_sqr_gens.append(sigma)
                elim.extend([sigma * tau for tau in elim])

        # Computing the real beta
        self._beta_first = dict()
        V = VectorSpace(GF(2), len(mod_sqr_gens))
        for v in V:
            for tau in sqr_gens:
                alpha_s = c_cor(tau, tau)
                s = tau^2
                for i in range(len(v)):
                    if v[i] == 1:
                        sigma_i = mod_sqr_gens[i]
                        alpha_s *= c_cor(sigma_i, s)
                        s = s * sigma_i
                self._beta_first[s] = alpha_s * beta_guess(s)

        # Correcting
        @cached_function
        def beta(sigma):
            s = galois_field_change(sigma, K)
            return self._beta_first[s]

        # Check:
        n = len(G)
        self._c_beta_error_index = list(G)
        self._c_beta_error = [0]*n
        c_beta = self._map_coboundary(beta)
        flag = False
        for i in range(n):
            sigma = self._c_beta_error_index[i]
            self._c_beta_error[i] = [0]*n
            for j in range(n):
                tau = self._c_beta_error_index[j]
                self._c_beta_error[i][j] = c_beta(sigma, tau) / c(sigma, tau)
                if self._c_beta_error[i][j] != 1:
                    flag = True
        if flag:
            raise Exception("The splitting field does not allow a splitting of c\n" + \
                            "See the function c_beta_error for more information.")

        return beta

    def c_beta_error(self, method, verbose=True, index=None, S=[]):
        if not self._is_cached('_c_beta_error'):
            try:
                self.splitting_map()
            except Exception:
                pass
        if self._is_cached('_beta') and 0 in self._beta:
            print "No error found!"
        if isinstance(method, str):
            if index is None:
                index = self._c_beta_error_index
                c_beta_error = self._c_beta_error
            else:
                conversion = [{k : v for v, k in enumerate(self._c_beta_error_index)}[s] for s in index]
                c_beta_error = [[self._c_beta_error[conversion[i]][conversion[j]] for j in range(len(index))]
                                                                                  for i in range(len(index))]
            index_dict = {k : v for v, k in enumerate(index)}
            n = len(index)
            method = method.lower()
            if method == 'table':
                if verbose:
                    print "Index set"
                    print list(enumerate(index))
                    print ""
                    l = max(floor(log(n-1,10))+2,3)
                    row_format=("{:>"+str(l)+"}")*(n+1)
                    print row_format.format("", *range(n))
                    for i in range(n):
                        print row_format.format(i, *[c_beta_error[i][j] for j in range(n)])
                return index, c_beta_error
            if method == 'relations':
                result = []
                for i,j in index_loop([n,n]):
                    ij = index_dict[index[i]*index[j]]
                    s = "a" + str(i) + " " + str(i) + "(a" + str(j) + ") = " +\
                        str(c_beta_error[i][j]) + " a" + str(ij)
                    if verbose:
                        print s
                    result.append(s)
                return result
            if method == 'congruences':
                K = index[0].parent().number_field()
                S1 = []
                for p in S:
                    S1.extend(K.primes_above(p))
                U = K.S_unit_group(S=S1)
                units = U.gens()
                orders = [u.order() for u in units]
                t = min([i for i in range(len(units)) if orders[i] in ZZ])
                N = orders[t]
                units = [K(u) for u in units]
                G_action = [matrix([U(index[i](u)).list() for u in units]) for i in range(n)]
                names = [0] * (n * len(units))
                for i in range(n):
                    for j in range(len(units)):
                        names[len(units)*i + j] = "a" + str(i) + "b" + str(j)
                big_mat = []
                big_tor_mat = []
                big_tor_val = []
                big_val = []
                for i,j in index_loop([n,n]):
                    ij = index_dict[index[i]*index[j]]
                    m = [[0]*len(names) for k in range(len(units))]
                    v = U(K(c_beta_error[i][j])).list()
                    for k in range(len(units)):
                        m[k][len(units)*i + k] += 1
                        m[k][len(units)*ij + k] -= 1
                        for l in range(len(units)):
                            m[l][len(units)*j + k] += G_action[i][k][l]
                    if verbose:
                        for k in range(len(units)):
                            s = ""
                            for l in range(len(names)):
                                if m[k][l] != 0:
                                    if len(s) > 0:
                                        s += " + "
                                    s += str(m[k][l]) + " " + names[l]
                            s += " = " + str(v[k])
                            if k == t:
                                s += " mod " + str(N)
                            print s
                    big_tor_mat.append(m[t])
                    del m[t]
                    big_mat.extend(m)
                    big_tor_val.append(v[t])
                    del v[t]
                    big_val.extend(v)
                return matrix(ZZ, big_mat), vector(ZZ, big_val), \
                       matrix(ZZ.quotient(N), big_tor_mat), \
                       vector(ZZ.quotient(N), big_tor_val), N, \
                       units
        print "No method", method, "found."
        print "Try 'table', 'relations' or 'congruences'"
    
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
    
    def splitting_map(self, index=0):
        if not self._is_cached('_beta') or 0 not in self._beta:
            self._beta = dict()
            self._beta[0] = self._first_splitting_map()
        if hasattr(index, "__iter__"):
            return [self.splitting_map(i) for i in index]
        if index in ZZ:
            if index not in self._beta:
                self._beta[index] = self._indexed_splitting_map(index)
            return self._beta[index]
        if index == 'all':
            return [self.splitting_map(i) for i in range(self.number_of_splitting_maps())]
        if index == 'conjugacy':
            return [self.splitting_map(i[0]) for i in self._conjugacy_determination()]
        raise Exception("Invalid index %s"%index)

    @cached_method
    def complete_field(self):
        self.fill_isogenies()
        Kl = self._Kl
        Kbeta = self.splitting_field()
        return composite_field(Kl, Kbeta)

    def _Kl_roots(self):
        self.fill_isogenies()
        Kl = self._Kl
        products = [1]
        result = []
        for tmp in Kl.subfields(degree=2):
            c = tmp[0].discriminant().squarefree_part()
            if c not in products:
                result.append(c)
                products.extend([(c*b).squarefree_part() for b in products])
        return result

    def cyclotomic_order(self):
        """Returns the smallest N such that Q(\zeta_N) contains the field of complete definition"""
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
        if not self._is_cached('_chi'):
            self._init_twist_characters()
        if hasattr(i, "__iter__"):
            return [self._twist_character_data(ii, j) for ii in i]
        if i in ZZ:
            if j == 1 and len(self._chi[i]) < 2:
                self._chi[i].append(dirichlet_to_galois(self._chi[i][0]))
            return self._chi[i][j]
        if i == 'all':
            return [self._twist_character_data(ii, j) for ii in range(self.number_of_splitting_maps())]
        if i == 'conjugacy':
            return [self._twist_character_data(ii[0], j) for ii in self._conjugacy_determination()]
        raise Exception("Invalid index %s."%i)
    
    def twist_character(self, index=0, galois=False):
        if galois:
            return self._twist_character_data(index, 1)
        else:
            return self._twist_character_data(index, 0)

    @cached_method
    def number_of_splitting_maps(self, count_conjugates=True):
        if count_conjugates:
            if not self._is_cached('_chi'):
                self._init_twist_characters()
            return len(self._chi)
        else:
            return len(self._conjugacy_determination)

    @cached_method
    def _conjugacy_determination(self):
        beta_ls = self.splitting_map("all")
        beta_del = [beta for beta in beta_ls]
        beta_dict = dict()
        for i in range(len(beta_ls)):
            beta_dict[beta_ls[i]] = i
        Kcomp = self.complete_field()
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
        return result

    @cached_method
    def conductor_restriction_of_scalars(self):
        self.fill_isogenies()
        K = self._E0.base_ring()
        Kcomp = self.complete_field()
        to_Kcomp = K.hom([a.minpoly().change_ring(Kcomp).roots()[0][0] for a in K.gens()], Kcomp)
        Ecomp = self._E0.change_ring(to_Kcomp)
        N_Ecomp = Ecomp.conductor()
        self._N_Res = Kcomp.discriminant()^2 * N_Ecomp.absolute_norm()
        # Proposition 1 of Milne, On the arithmetic of Abelian varieties
        return self._N_Res

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

### Example 1 from Quer
t = 3 # Parameter 1/5, 2, 27/25, 16/17, 125/121, 80/81, 1024/1025, 3024/3025, 250000/250001 are CM
K = field_with_root(QQ, t, names='sqrt_t')
t = K(t)
sqrt_t = sqrt(t)
a4 = -3 * sqrt_t * (4 + 5*sqrt_t)
a6 = 2 * sqrt_t * (2 + 14*sqrt_t + 11*t)
E = EllipticCurve([a4,a6])
if E.has_cm():
    print " Warning: A CM elliptic curve has been chosen!!!"
print E

Q = Quer_invariants(E); Q.add_isogenies(3)

### Twists of example 1
twist_ls = compute_possible_twists(K, K.prime_above(2))
for i in range(len(twist_ls)):
    print i, twist_elliptic_curve(E, twist_ls[i]).conductor().factor()
    print i, twist_elliptic_curve(E, sqrt_t * twist_ls[i]).conductor().factor()

### My twist
Q1 = Quer_invariants(twist_elliptic_curve(E, sqrt_t)); Q1.add_isogenies(3)

### Finding a small level
Kbig = CyclotomicField(24)
to_Kbig = K.hom([a.minpoly().change_ring(Kbig).roots()[0][0] for a in K.gens()], Kbig)
Ebig = E.change_ring(to_Kbig)
Qbig = Quer_invariants(Ebig); Qbig.add_isogenies(1); Qbig.add_isogenies(3);

### Using magma to compute with newforms
eps = magma.DirichletGroup(2^4 * 3^5)(product(magma.DirichletGroup(12).gens()))
cfs = magma.CuspForms(eps)
nfs1 = magma.Newforms(cfs)
nfs = []
for tmp in nfs1:
    f = tmp[1]
    K1 = magma.BaseField(f).sage()
    K2 = Qbig.splitting_image_field()
    if K1.degree() == K2.degree() and K1.discriminant().squarefree_part() == K2.discriminant().squarefree_part():
        nfs.append(f)

### Comparing L-functions to find the right newform
L = Qbig.splitting_image_field()
R.<T> = L[]
chi_ls = Qbig.twist_character('conjugacy')
Kf_ls = Qbig.splitting_image_field('conjugacy')
index_ls = range(len(nfs))
for p in prime_range(5, 100):
    Lchi_E_p = get_Lchi_E_p(Ebig, p, T)
    remove_ls = []
    for i in index_ls:
        f = nfs[i]
        Lchi_f_p = get_Lchi_f_p(f, chi_ls, Kf_ls, p, T)
        if Lchi_f_p != Lchi_E_p:
            remove_ls.append(i)
    for i in remove_ls:
        index_ls.remove(i)
print [nfs[i] for i in index_ls]

### Sanity check
for p in prime_range(5, 10):
    Lchi_E_p = get_Lchi_E_p(Ebig, p, T)
    print Lchi_E_p
    for f in nfs:
        Lchi_f_p = get_Lchi_f_p(f, chi_ls, Kf_ls, p, T)
        if Lchi_f_p == Lchi_E_p:
            print Lchi_f_p, "*"
        else:
            print Lchi_f_p
    print ""

### Using the twists
E1 = twist_elliptic_curve(E, twist_ls[3])
Q1 = Quer_invariants(E1); Q1.add_isogenies(3)
E2 = twist_elliptic_curve(E, sqrt_t)
Q2 = Quer_invariants(E2); Q2.add_isogenies(3);

### The curve:
t = 2 # Parameter 3, 9/2, 11/3, 25/6, 75/19, 289/72, 675/169, 1089/272, 31211/7803 are CM
s = 1 + 2 * t
K = composite_field(field_with_root(QQ, t, names='sqrt_t'), field_with_root(QQ, s, names='sqrt_s'))
sqrt_t = sqrt(K(t))
sqrt_s = sqrt(K(s))
a4 = -6*s*t*(5 + 5*sqrt_s + 10*sqrt_t + 5*t + 2*sqrt_s*sqrt_t)
a6 = 8 * (sqrt_s*sqrt_t)^3 * (1 + sqrt_t) * (7 + 15*sqrt_s + 14*sqrt_t + 7*t + 6*sqrt_s*sqrt_t)
E_1 = EllipticCurve([a4,a6])
if E_1.has_cm():
    print "Warning: A CM elliptic curve has been chosen!!!"
print E_1

### The twists
twist_ls = compute_possible_twists(K, K.prime_above(2))
E_2 = twist_elliptic_curve(E_1, twist_ls[-1])
E_3 = twist_elliptic_curve(E_1, twist_ls[-2])
E_1 = twist_elliptic_curve(E_1, sqrt_t - t)

### Constructing corresponding classes
Q1 = Quer_invariants(E_1); Q1.add_isogenies(2); Q1.add_isogenies(3)
#Q2 = Quer_invariants(E_2); Q2.add_isogenies(2); Q2.add_isogenies(3)
#Q3 = Quer_invariants(E_3); Q3.add_isogenies(2); Q3.add_isogenies(3)

### L-function functions
def get_Lchi_E_P(E, P, T):
    NP = P.absolute_norm()
    eP = NP.ord(p)
    if E.has_good_reduction(P):
        nP = E.reduction(P).count_points()
        aP = NP + 1 - nP
        bP = NP
    elif E.has_split_multiplicative_reduction(P):
        aP = 1
        bP = 0
    elif E.has_nonsplit_multiplicative_reduction(P):
        aP = -1
        bP = 0
    elif E.has_additive_reduction(P):
        aP = 0
        bP = 0
    else:
        raise Exception("Curve %s has no type of reduction at %s"%(E, P))
    return T^(2*eP) - aP * T^eP + bP

def get_Lchi_E_p(E, p, T):
    return product([get_Lchi_E_P(E, P, T) for P in E.base().primes_above(p)])

def get_Lchi_p_factors(E, p, L):
    R.<T> = L[]
    factors = dict()
    for (f, e) in get_Lchi_E_p(E, p, T).factor():
        factors[f] = e
    return factors

def sigma_to_poly(sigma, f):
    return f.parent()([sigma(c) for c in f.coefficients(sparse=False)])

def get_Lchi_p_invariant_factors(E, p, L):
    result = []
    factors = get_Lchi_p_factors(E, p, L)
    while len(factors) > 0:
        f,e = factors.items()[0]
        result.append(f)
        for sigma in L.galois_group():
            sigma_f = sigma_to_poly(sigma, f)
            factors[sigma_f] -= 1
            if factors[sigma_f] <= 0:
                factors.pop(sigma_f)
    return result

### L-function functions for magma newforms
def get_Lchi_f_p_chi(ap_f, eps_p, chi_p, Kf, p, T):
    """Gives the p-factor corresponding to chi f"""
    K1, ap_map, chi_map = composite_field(ap_f.parent(), chi_p.parent(), give_maps=True)
    a = ap_map(ap_f) * chi_map(chi_p)
    K2, chi_map, eps_map = composite_field(chi_p.parent(), eps_p.parent(), give_maps=True)
    b = p * chi_map(chi_p)^2 * eps_map(eps_p)
    a = a.minpoly().change_ring(Kf).roots()[0][0]
    b = a.minpoly().change_ring(Kf).roots()[0][0]
    L = T.parent().base_ring()
    result = T.parent().one()
    for (c,e) in Kf.gen().minpoly().change_ring(L).roots():
        sigma = Kf.hom([c], L)
        result *= T^2 - sigma(a)*T + sigma(b)
    return result

def get_Lchi_f_p(f, chi_ls, Kf_ls, p, T):
    """Gives the p-part of the polynomial defining the L-function corresponding to the newform f"""
    if f.parent() == magma:
        ap_f = magma.Coefficient(f, p).sage()
        eps_p = (magma.DirichletCharacter(f)(p)).sage()
        if eps_p in ZZ:
            eps_p = QQ(eps_p)
    else:
        ap_f = f.coefficient(p)
        eps_p = f.character()(p)
    result = T.parent().one()
    for i in range(len(chi_ls)):
        chi = chi_ls[i]
        Kf = Kf_ls[i]
        result *= get_Lchi_f_p_chi(ap_f, eps_p, chi(p), Kf, p, T)
    return result
        
### Analyzing the factors
# for p in prime_range(6,30):
#     print p
#     print get_Lchi_p_factors(Ecomp, p, composite_field(Lbeta, CyclotomicField(2*eps.order())))
#     print get_Lchi_p_invariant_factors(Ecomp, p, composite_field(Lbeta, CyclotomicField(2*eps.order())))

### Beta c check:
def beta_c_check(Q):
    Kcomp = Q.complete_field()
    G = Kcomp.galois_group()
    beta = Q.splitting_map()
    c = Q.c
    for sigma in G:
        print '   c:', [c(sigma, tau) for tau in G]
        print 'beta:', [beta(sigma) * beta(tau) * beta(sigma*tau)^(-1) for tau in G]
        print ''

### Making a usefull list
Kcomp = Q1.complete_field()
G = Kcomp.galois_group()
beta = Q1.splitting_map()
c = Q1.c
gens = G.gens()
orders = [g.order() for g in gens]
g2_ls = [product([gens[k]^j[k] for k in range(len(gens))]) for j in index_loop(orders)]
for i in index_loop(orders):
    g1 = product([gens[k]^i[k] for k in range(len(gens))])
    beta_ls = [beta(g1) * beta(g2) * beta(g1*g2)^(-1) for g2 in g2_ls]
    c_ls = [c(g1, g2) for g2 in g2_ls]
    print [beta_ls[k] / c_ls[k] for k in range(len(g2_ls))]
    print [b(g1) * g1(b(g2)) / b(g1*g2) for g2 in g2_ls]
    print ""

### tmp
def b(s):
    if s == sigma or s == sigma^0:
        return 1
    elif s == tau:
        return alpha
    else:
        return -alpha

### for real
def b(s):
    if s == sigma^0:
        return sqrt( s(gamma)/gamma )
    else:
        return -sqrt( s(gamma)/gamma )
    
### Trying to find the error
K = Q1.splitting_field()
G = K.galois_group()
beta = Q1.splitting_map()
c = Q1.c
ls = []
for a in range(2):
    for b in range(4):
        ls.append((a,b))
G = [product([G.gens()[i]^ab[i] for i in range(len(ab))]) for ab in ls]
print "      ", ls
for i in range(len(ls)):
    sigma = G[i]
    nls = [QQ(c(sigma, tau) / (beta(sigma) * beta(tau) * beta(sigma*tau)^(-1))) for tau in G]
    prls = []
    for n in nls:
        if n == -1:
            prls.append(" -1 ")
        else:
            prls.append("  1 ")
    print ls[i], prls

### tmp
Kbig = Q.complete_field()
G = Kbig.galois_group()
beta_f = Q._splitting_map_first_guess()
c = Q.c

def beta(s):
    if s in [index0[1],index0[5]]:
        return -beta_f(s)
    else:
        return beta_f(s)

# Check:
n = len(G)
Q._c_beta_error_index = list(G)
Q._c_beta_error = [0]*n
c_beta = Q._map_coboundary(beta)
flag = False
for i in range(n):
    sigma = Q._c_beta_error_index[i]
    Q._c_beta_error[i] = [0]*n
    for j in range(n):
        tau = Q._c_beta_error_index[j]
        Q._c_beta_error[i][j] = c_beta(sigma, tau) / c(sigma, tau)
        if Q._c_beta_error[i][j] != 1:
            flag = True

