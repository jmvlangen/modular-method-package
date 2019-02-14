from sage.combinat.combination import Combinations

class power_analyzer(SageObject):
    r"""
    A class that is useful in the analyzation of
    equations of the form

    MATH::

     f(x_1, \ldots x_n) = y^l

    where $f$ is a polynomial, $l$ is a strictly
    positive integer and $x_1, \ldots, x_n, y$ are
    variables which take on integer values.
    
    Throughout the documentation of this class we
    will call $f$ the polynomial.
    """

    def __init__(self, f):
        r"""
        INPUT:

        - ``f`` -- The polynomial, which should
          be defined over the rationals or a
          number field.
        """
        self._f = f

    def factor(self, K):
        r"""
        Gives a factor of the polynomial over a
        field extension K.

        INPUT:

        - ``K`` -- An extension of the field over
          which the polynomial is defined.

        OUTPUT:

        A factor of the polynomial  over the
        given field K.
        """
        return self._f.change_ring(K).factor()[0][0]

    def coprimality(self, K, g=None):
        r"""
        Determines which factors of the polynomial are
        coprime.

        INPUT:

        - ``K`` -- An extension of the field over which
          this equation is defined that is galois.
        - ``g`` -- A factor of the polynomial defined
          over the field K. By default will be the
          outcome of :meth:`factor`

        OUTPUT:
        
        A dictionary indexed by elements of the galois
        group of K, such that the s-th entry of the
        dictionary contains a number divisible by all
        primes that divide both g and the conjugate
        of g by s.
        """
        if g is None:
            g = self.factor(K)
        if g.is_homogeneous():
            return {s : g.macaulay_resultant(g.change_ring(s.as_hom())) for s in K.galois_group()}
        else:
            return {s : g.resultant(g.change_ring(s.as_hom())) for s in K.galois_group()}

    def relations(self, K, g=None):
        r"""
        Gives possible relations between the polynomial
        and a factor.
        
        Since g is a factor of the polynomial over some
        galois field extension, multiplying certain
        galois conjugates of g will result in a constant
        times the polynomial. A relation is a combination
        of the required galois homomorphism and the
        corresponding coefficient.

        INPUT:

        - ``K`` -- An extension of the field over which
          this equation is defined that is galois.
        - ``g`` -- A factor of the polynomial defined
          over the field K. By default will be the
          outcome of :meth:`factor`

        OUTPUT:

        A list of tuples of the form (s1, ..., sn, c)
        with s1, ..., sn elements of the galois group
        of K and c an element of K, such that the
        product of g conjugated by each si is equal
        to c times the polynomial.
        """
        if g is None:
            g = self.factor(K)
        n = ZZ(self._f.degree() / g.degree())
        G = K.galois_group()
        G_ls = [s for s in G if s != G(1)]
        result = []
        for s_ls in Combinations(G_ls, n - 1):
            candidate = g * product(g.change_ring(s.as_hom()) for s in s_ls)
            test = (candidate / self._f)
            if test.denominator() == 1:
                # We found a relation
                result.append(([G(1)] + list(s_ls), K(test.numerator().constant_coefficient())))
        return result

    @cached_method
    def bad_primes(self, K, g=None):
        r"""
        Gives the primes at which the factor g might not
        be an l-th power.

        INPUT:

        - ``K`` -- An extension of the field over which
          this equation is defined that is galois.
        - ``g`` -- A factor of the polynomial defined
          over the field K. By default will be the
          outcome of :meth:`factor`

        OUTPUT:

        A list of all prime ideals such that the
        valuation of g at those primes is not
        necessarily a multiple of l. Here we assume
        that g is evaluated at values such that the
        polynomial at those values is an l-th power.
        """
        if g is None:
            g = self.factor(K)
        coprime = self.coprimality(K, g)
        result = reduce(lambda I, J: I + J,
                        [reduce(lambda I, J: I.intersection(J),
                                [K.ideal(coprime[t*s^(-1)])
                                 for s, t in Combinations(s_ls, 2)],
                                K.ideal(C))
                         for s_ls, C in self.relations(K, g)])
        if result == K.ideal(0):
            return 0
        else:
            return result.prime_factors()

    @cached_method
    def constant_coeff(self, K, g, bad_primes_val, l):
        r"""
        Gives the constant such that g is an l-th power
        times that constant.

        INPUT:

        - ``K`` -- An extension of the field over which
          this equation is defined that is galois.
        - ``g`` -- A factor of the polynomial defined
          over the field K.
        - ``bad_primes_val`` -- A tuple of integers,
          where each integer is the valuation of g
          at the corresponding prime in the list
          returned by :meth:`bad_primes`.
        - ``l`` -- A strictly positive integer.

        OUTPUT:

        An element c of K such that we have
        ..MATH::

          g = c * u * x^l

        for some x in K and some unit u in the ring of
        integers of K. For this we assume that g is
        evaluated at values for which the polynomial
        is an l-th power.
        """
        bad_primes = self.bad_primes(K, g)
        I = product(bad_primes[i]^bad_primes_val[i] for i in range(len(bad_primes)))
        C = K.class_group()
        I_exp = C(I).list()
        # We find the class CJ such that I = CJ^l
        CJ = C(1)
        for i in range(len(C.gens())):
            g, x, y = xgcd(l, C.gens()[i].order())
            CJ = CJ * C.gens()[i]^(ZZ(I_exp[i] * x / g))
        J = (CJ^(-1)).ideal()
        return (I * (J^l)).gens_reduced()[0]
    
    def unit_relation(self, K, g, cg, relation, l):
        r"""
        Determines the relation on units defined by
        a relation on a factor of the polynomial.

        INPUT:

        - ``K`` -- A field extension of the field
          over which the polynomial is defined that
          is galois.
        - ``g`` -- A factor of the polynomial over
          the field K.
        - ``cg`` -- An element of K, such that
          g = cg * u * x^l for some x in K and some
          unit u of the ring of integers of K
        - ``relation`` -- A tuple of the form
          (s0, ..., sn, c) with s0, ..., sn elements
          of the galois group of K and c an element
          of K, such that the product of the
          conjugates of g by each si is equal to
          the polynomial times c.
        - ``l`` -- A prime number, the
          exponent to be considered whenever
          l-th powers are mentioned.

        OUTPUT:

        A matrix M and a vector v with integer
        coefficients, such that the vectors x
        satisfying M y = v are precisely the
        exponents y = (y1, ..., yn) such that
        the unit u = u1^y1 * ... * un^yn in
        g = cg * u * x^l satisfies the given
        relation on g. Here u1, ..., un are 
        the generators of the unit group of K
        as returned by :meth:`gens`.
        """
        s_ls, C = relation
        constant = C / product(s(cg) for s in s_ls)
        l_part = product(P^(ZZ(e/l)) for P, e in K.ideal(constant).factor()).gens_reduced()[0]
        U = K.unit_group()
        v = vector(U(constant / (l_part^l)).list()).column()
        galois_action = galois_on_units(K)
        M = sum(galois_action[s] for s in s_ls)
        relevant_units = [i for i in range(len(U.gens())) if (U.gens()[i].order() == Infinity or
                                                              l.divides(U.gens()[i].order()))]
        M = copy(M[relevant_units])
        v = copy(v[relevant_units])
        return M, v

    @cached_method
    def unit_coeffs(self, K, g, l, cg=None, bad_primes_val=None):
        r"""
        Gives the possible exponents of units u such
        that g = cg * u * x^l for a fixed constant
        cg in K and some x in K.

        INPUT:
        
        - ``K`` -- A field extension of the field
          over which the polynomial is defined that
          is galois.
        - ``g`` -- A factor of the polynomial over
          the field K.
        - ``l`` -- A prime number, the
          exponent to be considered whenever
          l-th powers are mentioned.
        - ``cg`` -- An element of K, such that
          g = cg * u * x^l for some x in K and some
          unit u of the ring of integers of K. By
          default will be set to the value returned
          by :meth:`constant_coeff`.
        - ``bad_primes_val`` -- A tuple of integers,
          where each integer is the valuation of g
          at the corresponding prime in the list
          returned by :meth:`bad_primes`. This
          argument only has to be provided if cg
          is not given.

        OUTPUT:
        
        A list of tuples of integers modulo l,
        consisting of all those tuples (x0, ..., xn)
        such that the unit u = u0^x0 * ... * un^xn 
        can satisfy g = cg * u * x^l for some x in K
        whilst agreeing with all relations on g given
        by :meth:`relations`. Here u0, ..., un are
        the generators of the unit group of K as
        returned by :meth:`gens`
        """
        if cg is None:
            if bad_primes_val is None:
                raise ValueError("Must provide the argument 'bad_primes_val' if the argument 'cg' is not provided.")
            cg = self.constant_coeff(K, g, bad_primes_val, l)
        relations = self.relations(K, g)
        n = len(relations)
        lin_rel = [self.unit_relation(K, g, cg, rel, l) for rel in relations]
        M, v = [block_matrix(n, 1, Mi) for Mi in zip(*lin_rel)]
        sol0 = M.change_ring(Integers(l)).solve_right(v)
        ker = M.change_ring(Integers(l)).right_kernel()
        V = ker.ambient_module()
        sol0 = V(vector(sol0))
        result = [sol0 + V(vi) for vi in ker]
        for vi in result:
            vi.set_immutable()
        return result

    def units(self, K, g, l, cg=None, bad_primes_val=None):
        r"""
        Gives the possible units u such that
        g = cg * u * x^l for a fixed constant
        cg in K and some x in K.

        INPUT:
        
        - ``K`` -- A field extension of the field
          over which the polynomial is defined that
          is galois.
        - ``g`` -- A factor of the polynomial over
          the field K.
        - ``l`` -- A prime number, the
          exponent to be considered whenever
          l-th powers are mentioned.
        - ``cg`` -- An element of K, such that
          g = cg * u * x^l for some x in K and some
          unit u of the ring of integers of K. By
          default will be set to the value returned
          by :meth:`constant_coeff`.
        - ``bad_primes_val`` -- A tuple of integers,
          where each integer is the valuation of g
          at the corresponding prime in the list
          returned by :meth:`bad_primes`. This
          argument only has to be provided if cg
          is not given.

        OUTPUT:
        
        A list of units u of the ring of integers
        of K that can satisfy g = cg * u * x^l for
        some x in K whilst agreeing with all
        relations on g given by :meth:`relations`.
        This list contains representatives for each
        such units modulo multiplication by l-th powers.
        """
        coeffs = self.unit_coeffs(K, g, l, cg=cg, bad_primes_val=bad_primes_val)
        U = K.unit_group()
        return [product(U.gens()[i]^ZZ(cf[i]) for i in range(len(U.gens()))) for cf in coeffs]

    def prime_data(self, K, g, l, primes, cg=None, bad_primes_val=None):
        r"""
        Gives the possible values of the variables of
        g modulo a prime p that can satisfy the
        relation g = cg * u * x^l for some x in K
        and some unit u.        

        INPUT:
        
        - ``K`` -- A field extension of the field
          over which the polynomial is defined that
          is galois.
        - ``g`` -- A factor of the polynomial over
          the field K.
        - ``l`` -- A prime number, the
          exponent to be considered whenever
          l-th powers are mentioned.
        - ``primes`` -- A list of prime numbers for
          which the possible values modulo that
          prime should be determined. It can also
          be set to a non-negative integer in which
          case this list will be initialized as all
          the prime numbers smaller than that number.
        - ``cg`` -- An element of K, such that
          g = cg * u * x^l for some x in K and some
          unit u of the ring of integers of K. By
          default will be set to the value returned
          by :meth:`constant_coeff`.
        - ``bad_primes_val`` -- A tuple of integers,
          where each integer is the valuation of g
          at the corresponding prime in the list
          returned by :meth:`bad_primes`. This
          argument only has to be provided if cg
          is not given.

        OUTPUT:
        
        A dictionary indexed by the primes in the
        given list of primes, such that the number of
        elements of a residue field of a prime above p
        in K is 1 modulo l. For each prime p the
        corresponding value is a dictionary indexed
        by possible units u such that g = cg * u * x^l
        for some x in K. For each such u the
        corresponding value is a list of tuples
        containing the possible values of the variables
        of g for which the above relation is satisfied
        modulo all primes above p in K.
        """
        if primes in ZZ:
            primes = prime_range(primes)
        if cg is None:
            if bad_primes_val is None:
                raise ValueError("Must provide the argument 'bad_primes_val' if the argument 'cg' is not provided.")
            cg = self.constant_coeff(K, g, bad_primes_val, l)
        u_dict = {cf : {} for cf in self.unit_coeffs(K, g, l, cg=cg)}
        U = K.unit_group()
        G = K.galois_group()
        n = len(g.variables())
        for p in primes:
            P = K.prime_above(p)
            if P in self.bad_primes(K, g):
                continue
            Fp = P.residue_field()
            if mod(len(Fp), l) != 1:
                continue
            xp = Fp.primitive_element()
            B = diagonal_matrix([Fp(u).log(xp) for u in U.gens()])
            done_cases = {}
            removed = []
            for cf in u_dict:
                case = B*cf
                case.set_immutable()
                if case in done_cases:
                    done_cf = done_cases[case]
                    if done_cf in removed:
                        removed.append(cf)
                    else:
                        u_dict[cf][p] = u_dict[done_cases[case]][p]
                else:
                    u = product(U.gens()[i]^ZZ(cf[i]) for i in range(len(cf)))
                    xgl_poly = [(g*(cg*u)^(-1)).change_ring(s.as_hom()).change_ring(Fp) for s in G]
                    arg_ls = tuple(tuple(arg) for arg in mrange([p]*n)
                                   if gcd(arg) != 0 and
                                   all(xgl_arg == 0 or l.divides(xgl_arg.log(xp))
                                       for xgl_arg in (xgl(arg) for xgl in xgl_poly)))
                    if len(arg_ls) == 0:
                        removed.append(cf)
                    else:
                        u_dict[cf][p] = arg_ls
                    done_cases[case] = cf
            for cf in removed:
                del u_dict[cf]
        if len(u_dict) == 0:
            return {}
        primes = u_dict[u_dict.keys()[0]].keys()
        return {p : {product(U.gens()[i]^ZZ(cf[i]) for i in range(len(cf))) :
                     u_dict[cf][p] for cf in u_dict}
                for p in primes}

    def prime_conditions(self, K, g, l, primes, cg=None, bad_primes_val=None):
        r"""
        Gives the possible values of the variables of
        g modulo a prime p that can satisfy the
        relation g = cg * u * x^l for some x in K
        and some unit u as a condition.        

        INPUT:
        
        - ``K`` -- A field extension of the field
          over which the polynomial is defined that
          is galois.
        - ``g`` -- A factor of the polynomial over
          the field K.
        - ``l`` -- A prime number, the
          exponent to be considered whenever
          l-th powers are mentioned.
        - ``primes`` -- A list of prime numbers for
          which the possible values modulo that
          prime should be determined. It can also
          be set to a non-negative integer in which
          case this list will be initialized as all
          the prime numbers smaller than that number.
        - ``cg`` -- An element of K, such that
          g = cg * u * x^l for some x in K and some
          unit u of the ring of integers of K. By
          default will be set to the value returned
          by :meth:`constant_coeff`.
        - ``bad_primes_val`` -- A tuple of integers,
          where each integer is the valuation of g
          at the corresponding prime in the list
          returned by :meth:`bad_primes`. This
          argument only has to be provided if cg
          is not given.

        OUTPUT:
        
        A dictionary indexed by the primes in the
        given list of primes, such that the number of
        elements of a residue field of a prime above p
        in K is 1 modulo l. For each prime p the
        corresponding value is a Condition such that
        the values of the variables satisfy this
        condition if and only if in that case
        g = cg * u * x^l modulo all primes above p in
        K for some unit u and some x in K. Note that
        the possible units u can be limited by
        their options modulo other primes.
        """
        data = self.prime_data(K, g, l, primes, cg=cg,
                               bad_primes_val=bad_primes_val)
        result = {}
        for p in data:
            T = pAdicTree(variables=self._f.variables(), prime=p, ring=QQ, full=False)
            Tr = T._root
            for val in {val for u in data[p] for val in data[p][u]}:
                Tr.children.add(pAdicNode(pAdics=Tr.pAdics(),
                                          coefficients=val,
                                          full=True))
            result[p] = TreeCondition(T)
        return result

def polynomial_split_on_basis(f, B):
    r"""
    Determines the part of the polynomial associated
    to each basis element.

    Given a (multivariate) polynomial over a
    number field K, one can see this as the sum of
    (multivariate) polynomials over the rationals
    times an associated basis element of a basis of
    K over the rationals. This method constructs
    these polynomials for a given basis.

    INPUT:

    - ``f`` -- A (multivariate) polynomial over some
      number field.
    - ``B`` -- A basis over QQ for the number field
      over which the polynomial is defined.

    OUTPUT:

    A list ls of polynomials over the rationals such
    that sum(ls[i] * B[i] for i in range(len(B))) == f.
    """
    K = f.parent().base_ring()
    M = matrix([K(B[i]).list() for i in range(len(B))]).transpose()
    if not M.is_invertible():
        raise ValueError("%s is not a basis for %s"%(B, K))
    R = f.parent().change_ring(QQ)
    # For each monomial m
    #   M^(-1) * vector(f.monomial_coefficient(m).list())
    # is the vector that expresses the corresponding coefficient
    # of f in terms of the basis B
    monomials = [[cf * R(m)
               for cf in (M^(-1) * vector(f.monomial_coefficient(m).list()))]
              for m in f.monomials()]
    return [sum(ls) for ls in zip(*monomials)]
