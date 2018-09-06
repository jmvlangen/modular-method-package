from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing

from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class

class FreyCurve(EllipticCurve_generic):
    r"""
    A Frey-curve.

    A Frey-curve is a curve defined by the (unknown) solution
    of some Diophantine equation. This is represented by an
    elliptic curve over some polynomial ring, where the
    coefficient ring is the actual field of definition of the
    curve and the variables represent the values of the
    solution.
    """
    def __init__(self, curve, solution_ring=ZZ, conversion=None, condition=None):
        r"""
        Constructor of a Frey curve

        INPUT:

        - ``curve`` -- An elliptic curve defined over some
          polynomial ring or any argument that would produce
          such a curve when passed to the constructor
          EllipticCurve. This will become the Frey curve.
        - ``solution_ring`` -- A ring (default: ZZ) that
          has a natural map into the base ring of the
          polynomial ring over which this elliptic curve
          is defined. This is the ring in which the variables
          of the polynomial ring over which this curve is
          defined can take values.
        - ``conversion`` -- A map (default: None) from the
          solution ring to the base ring of the polynomial
          ring over which the coefficients of this curve
          are defined. If set to None, will attempt to find
          such a map by trying in this order maps given by:
           - coerce_map_from
           - convert_map_from
           - Hom( , ).an_element()
        - ``condition`` -- A Condition object or None
          (default: None) giving a condition which must hold
          on the variables in this Frey-curve. If set to
          None will assume that all values for these variables
          are possible instead.
        """
        if not isinstance(curve, EllipticCurve_generic):
            curve = EllipticCurve(curve)
        S = curve.base_ring()
        if not is_PolynomialRing(S) or is_MpolynomialRing(S):
            raise ValueError("The coefficient ring %s is not a polynomial ring."%S)
        base = S.base_ring()
        EllipticCurve_generic.__init__(self, S, curve.a_invariants())
        self._R = solution_ring
        if conversion is None:
            conversion = base.coerce_map_from(solution_ring)
        if conversion is None:
            conversion = base.convert_map_from(solution_ring)
        if conversion is None:
            conversion = Hom(solution_ring, base).an_element()
        self._R_to_base = conversion
        self._condition = condition

    def definition_ring(self):
        r"""
        Gives the ring over which this Frey curve is defined.

        Even though the base ring of a Frey curve is a polynomial
        ring, it is assumed that the Frey curve is defined over
        the base ring of this polynomial ring, since the variables
        of this polynomial ring are assumed to take on unknown
        values in this ring.

        OUTPUT:
        
        The base ring of the polynomial ring over which this
        Frey curve is defined.
        """
        return self.base_ring().base()
    
    def parameters(self):
        r"""
        Gives the parameters on which this Frey curve depends.

        OUTPUT:

        The generators of the polynomial ring over which this
        Frey curve is defined.
        """
        return self.base().gens()

    @cached_method
    def primes_of_possible_additive_reduction(self):
        r"""
        Computes the primes at which this curve could have additive reduction.

        OUTPUT:

        A list of primes, i.e. prime numbers if the base ring of this
        frey curve is a polynomial ring over a subring of QQ or maximal
        ideals if it is a polynomial ring over a subring of a number
        field. This list aims to contain all primes at which this Frey
        curve might have additive reduction, but could not be sufficient.
        In case of uncertainty, will print a warning.
        """
        R = self.base()
        n = len(R.gens())
        K = R.base()
        c4 = self.c4()
        D = self.discriminant()
        if K.is_subring(QQ):
            if n == 1:
                return QQ(c4.resultant(D)).numerator().prime_factors()
            elif n == 2:
                return QQ(c4.macaulay_resultant(D)).numerator().prime_factors()
            else:
                N = lcm(lcm(QQ(c).denominator() for c in c4.coefficients()),
                        lcm(QQ(c).denominator() for c in D.coefficients()))
                M = gcd(gcd(ZZ(c) for c in (N * c4).coefficients()),
                        gcd(ZZ(c) for c in (N * D).coefficients()))
                return QQ(M/N).numerator().prime_factors()
        if n == 1:
            return K.ideal(c4.resultant(D)).prime_factors()
        elif n == 2 and c4.is_homogeneous() and D.is_homogeneous():
            return K.ideal(c4.macaulay_resultant(D)).prime_factors()
        else:
            I = sum(K.ideal(c) for c in c4.coefficients())
            J = sum(K.ideal(c) for c in D.coefficients())
            result = (I + J).prime_factors()
            print "Warning: Assuming that %s and %s"%(c4,D) + \
                  "are coprime outside %s."%(tuple(P._repr_short() for P in result),)
            return result

    @cached_method
    def _initial_tree(self, prime, verbose=False):
        r"""
        Gives the tree of possible values of the parameters for
        a given prime.

        INPUT:

        - ``prime`` -- A (maximal) prime ideal of the ring in
          which the parameters take value or any generator
          thereof if it is principal.
        - ``verbose`` -- A boolean (default: False) indicating
          whether additional information should be printed
          during computation.
        
        OUTPUT:

        A pAdicTree containing all the possible values of the
        the parameters in the completion at the given prime.
        """
        pAdics = pAdicBase(self.definition_ring(), prime)
        Tfull = pAdicTree(variables=self.parameters(),
                          prime=pAdics.prime_below(self._R))
        return self._condition.pAdic_tree(pAdic_tree=Tfull, verbose=verbose)

    @cached_method(key=lambda self, prime, verbose, precision_cap: (prime, precision_cap))
    def local_data(self, prime, verbose=False, precision_cap=20):
        r"""
        Gives a minimal model of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime, i.e. a prime number if this Frey
          curve is defined over a subring of QQ or a maximal ideal
          of the ring of integers if it is defined over a subring
          of a number field.
        - ``verbose`` -- A boolean (default: False) indicating
          whether additional information should be printed
          during computation.
        - ``precision_cap`` A strictly positive integer
          (default: 20) giving the maximal precision level to be
          used in p-Adic arithmetic.

        OUTPUT:
        
        An elliptic curve that is a minimal model of this curve
        at the given prime. This could be a ConditionalValue as
        the minimal model might depend on the value of the
        parameters in this Frey curve.
        """
        pAdics = pAdicBase(self.definition_ring(), prime)
        Tp = _initial_tree(pAdics.prime_below(self._R),
                           verbose=verbose)
        return performTatesAlgorithm(self,
                                     initial_values=Tp,
                                     coefficient_ring=self.base(),
                                     pAdics=pAdics,
                                     verbose=verbose,
                                     precision_cap=precision_cap)

    @cached_method(key=lambda self, prime, verbose, precision_cap: (prime, precision_cap))
    def minimal_model(self, prime, verbose=False, precision_cap=20):
        r"""
        Gives a minimal model of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime, i.e. a prime number if this Frey
          curve is defined over a subring of QQ or a maximal ideal
          of the ring of integers if it is defined over a subring
          of a number field.
        - ``verbose`` -- A boolean (default: False) indicating
          whether additional information should be printed
          during computation.
        - ``precision_cap`` A strictly positive integer
          (default: 20) giving the maximal precision level to be
          used in p-Adic arithmetic.

        OUTPUT:
        
        An elliptic curve that is a minimal model of this curve
        at the given prime. This could be a ConditionalValue as
        the minimal model might depend on the value of the
        parameters in this Frey curve.
        """
        if self.local_data.is_in_cache(prime, verbose, precision_cap):
            local_data = self.local_data(prime, verbose, precision_cap)
            if isinstance(local_data, FreyCurveLocalData):
                return local_data.minimal_model()
            result = []
            for data, T in local_data:
                flag = True
                for i in range(len(result)):
                    if result[i][0] == data.minimal_model():
                        result[i] = (data.minimal_model(),
                                     result[i][1].union(T))
                        flag = False
                        break
                if flag:
                    result.append((data.minimal_model(), T))
            return ConditionalValue(result)
            
        pAdics = pAdicBase(self.definition_ring(), prime)
        Tp = _initial_tree(pAdics.prime_below(self._R),
                           verbose=verbose)
        return performTatesAlgorithm(self,
                                     initial_values=Tp,
                                     coefficient_ring=self.base(),
                                     pAdics=pAdics,
                                     verbose=verbose,
                                     precision_cap=precision_cap,
                                     only_calculate=['minimal_model'])

    @cached_method(key=lambda self, prime, verbose, precision_cap: (prime, precision_cap))
    def kodaira_symbol(self, prime, verbose=False, precision_cap=20):
        r"""
        Gives the kodaira symbol of the reduction at a given prime.

        INPUT:

        - ``prime`` -- A prime, i.e. a prime number if this Frey
          curve is defined over a subring of QQ or a maximal ideal
          of the ring of integers if it is defined over a subring
          of a number field.
        - ``verbose`` -- A boolean (default: False) indicating
          whether additional information should be printed
          during computation.
        - ``precision_cap`` A strictly positive integer
          (default: 20) giving the maximal precision level to be
          used in p-Adic arithmetic.

        OUTPUT:
        
        The KodairaSymbol representing the type of reduction of
        this curve at the given prime. This could be a
        ConditionalValue as it might depend on the value of
        the parameters in this curve.
        """
        if self.local_data.is_in_cache(prime, verbose, precision_cap):
            local_data = self.local_data(prime, verbose, precision_cap)
            if isinstance(local_data, FreyCurveLocalData):
                return local_data.kodaira_symbol()
            result = []
            for data, T in local_data:
                flag = True
                for i in range(len(result)):
                    if result[i][0] == data.kodaira_symbol():
                        result[i] = (data.kodaira_symbol(),
                                     result[i][1].union(T))
                        flag = False
                        break
                if flag:
                    result.append((data.kodaira_symbol(), T))
            return ConditionalValue(result)
        
        pAdics = pAdicBase(self.definition_ring(), prime)
        Tp = _initial_tree(pAdics.prime_below(self._R),
                           verbose=verbose)
        return performTatesAlgorithm(self,
                                     initial_values=Tp,
                                     coefficient_ring=self.base(),
                                     pAdics=pAdics,
                                     verbose=verbose,
                                     precision_cap=precision_cap,
                                     only_calculate=['minimal_model'])
    
    @cached_method(key=lambda self, prime, verbose, precision_cap: (prime, precision_cap))
    def conductor_exponent(self, prime, verbose=False, precision_cap=20):
        r"""
        Gives the conductor exponent at a given prime.

        INPUT:

        - ``prime`` -- A prime, i.e. a prime number if this Frey
          curve is defined over a subring of QQ or a maximal ideal
          of the ring of integers if it is defined over a subring
          of a number field.
        - ``verbose`` -- A boolean (default: False) indicating
          whether additional information should be printed
          during computation.
        - ``precision_cap`` A strictly positive integer
          (default: 20) giving the maximal precision level to be
          used in p-Adic arithmetic.

        OUTPUT:

        The exponent of the conductor of this Frey curve at
        the given prime. This could be a ConditionalValue as
        the conductor exponent might depend on the value of
        the parameters in this Frey curve.
        """
        if self.local_data.is_in_cache(prime, verbose, precision_cap):
            local_data = self.local_data(prime, verbose, precision_cap)
            if isinstance(local_data, FreyCurveLocalData):
                return local_data.conductor_valuation()
            result = []
            for data, T in local_data:
                flag = True
                for i in range(len(result)):
                    if result[i][0] == data.conductor_valuation():
                        result[i] = (data.conductor_valuation(),
                                     result[i][1].union(T))
                        flag = False
                        break
                if flag:
                    result.append((data.conductor_valuation(), T))
            return ConditionalValue(result)
        
        pAdics = pAdicBase(self.definition_ring(), prime)
        Tp = _initial_tree(pAdics.prime_below(self._R),
                           verbose=verbose)
        return performTatesAlgorithm(self,
                                     initial_values=Tp,
                                     coefficient_ring=self.base(),
                                     pAdics=pAdics,
                                     verbose=verbose,
                                     precision_cap=precision_cap,
                                     only_calculate=['conductor'])

    def conductor(self, additive_primes=None, verbose=False,
                  precision_cap=20):
        r"""
        Computes the possible conductors of this Frey curve.

        INPUT:

        - ``additive_primes`` -- An iterable containing prime ideals
          or prime numbers, if the field of definition is QQ, that
          contains all the primes at which this curve can have
          additive reduction. If set to None will compute this
          by using the method primes_of_possible_additive_reduction
        - ``verbose`` -- A boolean (default: False) indicating
          whether additional information should be printed
          during computation. If set to True will print
          information about the current step of computation
          during the computation.
        - ``precision_cap`` A strictly positive integer
          (default: 20) giving the maximal precision level to be
          used in p-Adic arithmetic.

        OUTPUT:

        A conditional expression that gives the conductor of
        this Frey curve for each possible value of the
        parameters on which it depends.
        """
        if additive_primes is None:
            additive_primes = self.primes_of_possible_additive_reduction()
        result = product(P^self.conductor_exponent(P, verbose=verbose,
                                                 precision_cap=precision_cap)
                         for P in additive_primes)
        result *= "Rad_P( " + str(self.discriminant().factor()) + " )"
        return result, additive_primes
        
class FreyCurveLocalData(EllipticCurveLocalData):    
    def __init__(self, elliptic_curve, prime,
                 conductor_valuation,
                 discriminant_valuation,
                 kodaira_symbol,
                 tamagawa_number,
                 reduction_type):
        self._set_elliptic_curve(elliptic_curve)
        self._fp = conductor_valuation
        self._val_disc = discriminant_valuation
        if isinstance(kodaira_symbol, KodairaSymbol_class):
            self._KS = kodaira_symbol
        else:
            self._KS = KodairaSymbol(kodaira_symbol)
        self._cp = tamagawa_number
        self._reduction_type = reduction_type
        
    def _set_elliptic_curve(self, elliptic_curve):
        self._Emin = elliptic_curve
        self._Emin_reduced = elliptic_curve
        
    def same_local_model(self, other):
        return isinstance(other, EllipticCurveLocalData) and \
        self.prime() == other.prime() and \
        self.kodaira_symbol() == other.kodaira_symbol() and \
        self.conductor_valuation() == other.conductor_valuation() and \
        self.tamagawa_number() == other.tamagawa_number()
        
    def same_elliptic_data(self, other):
        return self.same_local_model(other) and \
        self.minimal_model == other.minimal_model()
    
    def __eq__(self, other):
        return isinstance(other, ParametrizedLocalData) and \
        self.same_elliptic_data(other)
        
    def __ne__(self, other):
        return not isinstance(other, ParametrizedLocalData) or \
        not self.same_elliptic_data(other)
