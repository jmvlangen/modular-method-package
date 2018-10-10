def get_newforms(level, character=None, algorithm='sage', minimal_coeffs=QQ):
    r"""
    Computes the newforms of a given level and character.

    INPUT:

    - ``level`` -- A strictly positive integer indicating
      the level of the newforms.
    - ``character`` -- A dirichlet character of which the
      conductor divides the given level.
    - ``algorithm`` -- One of the following possible
      arguments:
       - 'sage' - (default) If sage should be used to
         compute the newforms.
       - 'magma' - If magma should be used to compute
         the newforms.
    - ``minimal_coeffs`` -- A number field or the rationals
      (default: QQ) that should be contained inf the
      coefficient field of each newform computed.
    
    OUTPUT:
    
    A list of instances of Newform_wrapped that contains
    exactly one newform in each galois orbit of newforms
    that are invariant under \gamma_1(level) with as
    character the given character. Furthermore the
    coefficient field of each of these newforms extends
    the given field minimal_coeffs.

    Note that these are instances of Newform_wrapped to
    provide a uniform way of computing with the newforms
    after they are created that is independant of the
    algorithm used to compute them.
    """
    if algorithm == 'sage':
        if character is None:
            nfs = Newforms(level)
        else:
            eps = character.primitive_character().extend(level)
            nfs = Newforms(eps)
        result = [Newform_wrapped_sage(f) for f in nfs]
    elif algorithm == 'magma':
        if character is None:
            cfs = magma.CuspForms(level)
            nfs = magma.Newforms(cfs)
        else:
            eps = character.primitive_character()
            gens = eps.parent().unit_gens()
            Dm = magma.DirichletGroup(eps.modulus(), magma(eps.base_ring()))
            candidate = False
            for eps_m in Dm.elements():
                candidate = all(eps_m(n) == eps(n) for n in gens)
                if candidate:
                    break
            if candidate: # We found a matching character
                eps_m = magma.DirichletGroup(level, magma(eps.base_ring()))(eps_m)
                cfs = magma.CuspForms(eps_m)
                nfs = magma.Newforms(cfs)
            else:
                raise ValueError("There is no dirichlet character in magma matching %s"%(eps,))
        result = [Newform_wrapped_magma(orbit[1]) for orbit in nfs]
    else:
        raise ValueError("%s is not a valid algorithm."%(algorithm,))
    if minimal_coeffs == QQ:
        return result
    else:
        poly = minimal_coeffs.gen().minpoly()
        return [f for f in result if len(poly.change_ring(f.coefficient_field()).roots()) > 0]

class Newform_wrapped(SageObject):
    r"""
    The wrapped version of a newform.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    This class is a base class, but should not be used
    as an instance. It rather provides a template for all
    classes that inherit it.
    """

    def level(self):
        r"""
        Gives the level of the newform.
        
        OUTPUT:

        A non-negative integer describing the level of this
        newform.
        """
        raise NotImplementedError()

    def character(self):
        r"""
        Gives the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform
        as a primitive character.
        """
        raise NotImplementedError()
        
    def coefficient(self, n):
        r"""
        Give the n-th coefficient of this newform.

        INPUT:
        
        - ``n`` -- A non-negative integer.

        OUTPUT:
        
        The n-th coefficient of the q-expansion of
        this newform at infinity.
        """
        raise NotImplementedError()

    def coefficient_field(self):
        r"""
        Gives the field over which the coefficients of
        this newform are defined.

        OUTPUT:

        The field over which the coefficients of the
        q-expansion of this newform at infinity are
        defined.
        """
        raise NotImplementedError()

    def q_expansion(self, prec=20):
        """
        Gives the q-expansion of this newform.

        INPUT:
        
        - ``prec`` -- A non-negative integer (default: 20)
          giving a bound on the powers that should be present
          in this q-expansion.

        OUTPUT:

        The q-expansion of this newform at infinity given
        as a power series in q with coefficients in the
        coefficient field of this newform and capped at
        the prec.
        """
        R.<q> = self.coefficient_field()[[]]
        result = sum(self.coefficient(n) * q^n for n in range(prec))
        return result.add_bigoh(prec)

    @cached_method
    def _trace_power_formula(self, power):
        R.<x,y> = QQ[]
        return polynomial_to_symmetric(x^power + y^power)

    def trace_of_frobenius(self, prime, power=1):
        """
        Gives the trace of frobenius under the galois
        representation associated to this newform.
        
        INPUT:

        - ``prime`` -- A prime number indicating the
          frobenius element to be used.
        - ``power`` -- A non-negative number (default: 1)
          If set to any value larger than 1, will compute
          the trace of the frobenius element to the given
          power instead.

        OUTPUT:

        The trace of rho_f(frob_p^n), where rho_f is the
        mod-l or l-adic galois representation associated
        to this newform, frob_p is the frobenius element
        at prime, and n is the given argument power.

        Since the result does not depend on the choice
        of l, this result will be an element of the
        coefficient field of the newform. The only
        condition is that l and p must be coprime.

        Will give a ValueError if the prime divides the
        level of this newform, since in that case all
        mentioned galois representations are ramified.
        """
        if prime not in ZZ or not prime.is_prime():
            raise ValueError("%s is not a prime number."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        if power == 1:
            return self.coefficient(prime)
        T = self.trace_of_frobenius(prime)
        D = self.determinant_of_frobenius(prime)
        return self._trace_power_formula(power)(T, D)

    def determinant_of_frobenius(self, prime, power=1):
        """
        Gives the determinant of frobenius under the
        galois representation associated to this newform.
        
        INPUT:

        - ``prime`` -- A prime number indicating the
          frobenius element to be used.
        - ``power`` -- A non-negative number (default: 1)
          If set to any value larger than 1, will compute
          the trace of the frobenius element to the given
          power instead.

        OUTPUT:

        The determinant of rho_f(frob_p^n), where rho_f is the
        mod-l or l-adic galois representation associated
        to this newform, frob_p is the frobenius element
        at prime, and n is the given argument power.

        Since the result does not depend on the choice
        of l, this result will be an element of the
        coefficient field of the newform. The only
        condition is that l and p must be coprime.

        Will give a ValueError if the prime divides the
        level of this newform, since in that case all
        mentioned galois representations are ramified.
        """
        if prime not in ZZ or not prime.is_prime():
            raise ValueError("%s is not a prime number."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        D = self.character()(prime) * prime
        return D^power

    def characteristic_polynomial(self, prime, power=1):
        """
        Gives the characteristic polynomial of the frobenius
        element acting on this newform.
        
        INPUT:

        - ``prime`` -- A prime number indicating the
          frobenius element to be used.
        - ``power`` -- A non-negative number (default: 1)
          If set to any value larger than 1, will compute
          the trace of the frobenius element to the given
          power instead.

        OUTPUT:

        The characteristic polynomial of rho_f(frob_p^n),
        where rho_f is the mod-l or l-adic galois
        representation associated to this newform,
        frob_p is the frobenius element at prime, and n
        is the given argument power.

        Since the result does not depend on the choice
        of l, this result will be an element of the
        coefficient field of the newform. The only
        condition is that l and p must be coprime.

        Will give a ValueError if the prime divides the
        level of this newform, since in that case all
        mentioned galois representations are ramified.
        """
        T = self.trace_of_frobenius(prime, power=power)
        D = self.determinant_of_frobenius(prime, power=power)
        K, T_map, D_map = composite_field(T.parent(), D.parent())
        R.<x> = K[]
        return x^2 - T_map(T)*x + D_map(D)
    
    qexp = q_expansion

    def _repr_(self):
        """
        Gives a string representation of self
        """
        return str(self.q_expansion())

    def _latex_(self):
        """
        Gives a latex representation of self.
        """
        return latex(self.q_expansion())

class Newform_wrapped_sage(Newform_wrapped):
    r"""
    The wrapped version of a newform in sage.
    """

    def __init__(self, newform):
        self._f = newform

    def level(self):
        r"""
        Gives the level of the newform.
        
        OUTPUT:

        A non-negative integer describing the level of this
        newform.
        """
        self._f.level()

    def character(self):
        r"""
        Gives the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform
        as a primitive character.
        """
        self._f.character()
        
    def coefficient(self, n):
        r"""
        Give the n-th coefficient of this newform.

        INPUT:
        
        - ``n`` -- A non-negative integer.

        OUTPUT:
        
        The n-th coefficient of the q-expansion of
        this newform at infinity.
        """
        self._f.coefficient(n)

    def coefficient_field(self):
        r"""
        Gives the field over which the coefficients of
        this newform are defined.

        OUTPUT:

        The field over which the coefficients of the
        q-expansion of this newform at infinity are
        defined.
        """
        self._f.base_ring()

    def q_expansion(self, prec=20):
        """
        Gives the q-expansion of this newform.

        INPUT:
        
        - ``prec`` -- A non-negative integer (default: 20)
          giving a bound on the powers that should be present
          in this q-expansion.

        OUTPUT:

        The q-expansion of this newform at infinity given
        as a power series in q with coefficients in the
        coefficient field of this newform and capped at
        the prec.
        """
        self._f.q_expansion(prec=20)

    def _repr_(self):
        """
        Gives a string representation of self.
        """
        return str(self._f)

    def _latex_(self):
        """Gives a latex representation of self."""
        return latex(self._f)

    qexp = q_expansion

class Newform_wrapped_magma(Newform_wrapped):
    r"""
    The wrapped version of a newform in magma.
    """

    def __init__(self, newform):
        self._f = newform

    def level(self):
        r"""
        Gives the level of the newform.
        
        OUTPUT:

        A non-negative integer describing the level of this
        newform.
        """
        self._f.Level().sage()

    @cached_method
    def character(self):
        r"""
        Gives the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform
        as a primitive character.
        """
        eps_f = self._f.DirichletCharacter()
        N = eps_f.Modulus().sage()
        N0 = eps_f.Conductor().sage()
        L = eps_f.BaseRing().sage()
        gens = Integers(N).unit_gens()
        for eps in DirichletGroup(N0, base_ring=L):
            if all(magma(eps(n)) == eps_f(n) for n in gens):
                return eps
        raise ValueError("No sage character corresponds to %s"%(eps_f,))
        
    def coefficient(self, n):
        r"""
        Give the n-th coefficient of this newform.

        INPUT:
        
        - ``n`` -- A non-negative integer.

        OUTPUT:
        
        The n-th coefficient of the q-expansion of
        this newform at infinity.
        """
        self._f.Coefficient(n).sage()

    def coefficient_field(self):
        r"""
        Gives the field over which the coefficients of
        this newform are defined.

        OUTPUT:

        The field over which the coefficients of the
        q-expansion of this newform at infinity are
        defined.
        """
        self._f.BaseField().sage()

    def _repr_(self):
        """Gives a string representation of self."""
        return str(self._f)

    def _latex_(self):
        """Gives a latex representation of self."""
        return latex(self._f)
