from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing
from sage.rings.morphism import RingHomomorphism_from_base

from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class

from sage.rings.number_field.number_field import is_NumberField

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
    def __init__(self, curve, parameter_ring=ZZ, conversion=None, condition=None):
        r"""
        Constructor of a Frey curve

        INPUT:

        - ``curve`` -- An elliptic curve defined over some
          polynomial ring or any argument that would produce
          such a curve when passed to the constructor
          EllipticCurve. This will become the Frey curve.
        - ``parameter_ring`` -- A ring (default: ZZ) that
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
        if not (is_PolynomialRing(S) or is_MPolynomialRing(S)):
            raise ValueError("The coefficient ring %s is not a polynomial ring."%S)
        base = S.base_ring()
        EllipticCurve_generic.__init__(self, S, curve.a_invariants())
        self._R = parameter_ring
        if conversion is None:
            conversion = base.coerce_map_from(parameter_ring)
        if conversion is None:
            conversion = base.convert_map_from(parameter_ring)
        if conversion is None:
            conversion = Hom(parameter_ring, base).an_element()
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

    def parameter_ring(self):
        r"""
        Gives the ring in which the parameters can take values.

        OUTPUT:

        A ring in which the parameters of this Frey curve take
        values. 
        """
        return self._R

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
            elif n == 2 and c4.is_homogeneous() and D.is_homogeneous():
                return QQ(c4.macaulay_resultant(D)).numerator().prime_factors()
            else:
                N = lcm(lcm(QQ(c).denominator() for c in c4.coefficients()),
                        lcm(QQ(c).denominator() for c in D.coefficients()))
                M = gcd(gcd([ZZ(c) for c in (N * c4).coefficients()]),
                        gcd([ZZ(c) for c in (N * D).coefficients()]))
                result = QQ(M/N).numerator().prime_factors()
                print "Warning: Assuming that %s and %s "%(c4,D) + \
                      "are coprime outside %s."%(tuple(result),)
                return result
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
        - ``verbose`` -- A boolean value or an integer
          (default: False). When set to True or any value
          larger then zero will print comments to stdout
          about the computations being done whilst busy. If
          set to False or 0 will not print such comments.
          If set to any negative value will also prevent
          the printing of any warnings.
          If this method calls any method that accepts an
          argument verbose will pass this argument to it.
          If such a method fulfills a minor task within
          this method and the argument verbose was larger
          than 0, will instead pass 1 less than the given
          argument. This makes it so a higher value will
          print more details about the computation than a
          lower one.
        
        OUTPUT:

        A pAdicTree containing all the possible values of the
        the parameters in the completion at the given prime.
        """
        Tfull = pAdicTree(variables=self.parameters(),
                          pAdics=pAdicBase(self._R, prime))
        if self._condition is None:
            return Tfull
        else:
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
        - ``verbose`` -- A boolean value or an integer
          (default: False). When set to True or any value
          larger then zero will print comments to stdout
          about the computations being done whilst busy. If
          set to False or 0 will not print such comments.
          If set to any negative value will also prevent
          the printing of any warnings.
          If this method calls any method that accepts an
          argument verbose will pass this argument to it.
          If such a method fulfills a minor task within
          this method and the argument verbose was larger
          than 0, will instead pass 1 less than the given
          argument. This makes it so a higher value will
          print more details about the computation than a
          lower one.
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
                           verbose=(verbose-1 if verbose>0 else verbose))
        result = performTatesAlgorithm(self,
                                       initial_values=Tp,
                                       coefficient_ring=self.base(),
                                       pAdics=pAdics,
                                       verbose=verbose,
                                       precision_cap=precision_cap)
        if len(result) == 1:
            return result[0][0]
        else:
            return result

    @cached_method(key=lambda self, prime, verbose, precision_cap: (prime, precision_cap))
    def minimal_model(self, prime, verbose=False, precision_cap=20):
        r"""
        Gives a minimal model of this curve at a given prime.

        INPUT:

        - ``prime`` -- A prime, i.e. a prime number if this Frey
          curve is defined over a subring of QQ or a maximal ideal
          of the ring of integers if it is defined over a subring
          of a number field.
        - ``verbose`` -- A boolean value or an integer
          (default: False). When set to True or any value
          larger then zero will print comments to stdout
          about the computations being done whilst busy. If
          set to False or 0 will not print such comments.
          If set to any negative value will also prevent
          the printing of any warnings.
          If this method calls any method that accepts an
          argument verbose will pass this argument to it.
          If such a method fulfills a minor task within
          this method and the argument verbose was larger
          than 0, will instead pass 1 less than the given
          argument. This makes it so a higher value will
          print more details about the computation than a
          lower one.
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
                           verbose=(verbose-1 if verbose>0 else verbose))
        result = performTatesAlgorithm(self,
                                       initial_values=Tp,
                                       coefficient_ring=self.base(),
                                       pAdics=pAdics,
                                       verbose=verbose,
                                       precision_cap=precision_cap,
                                       only_calculate=['minimal_model'])
        if len(result) == 1:
            return result[0][0][0]
        else:
            return ConditionalValue([(val[0], con) for val, con in result])

    @cached_method(key=lambda self, prime, verbose, precision_cap: (prime, precision_cap))
    def kodaira_symbol(self, prime, verbose=False, precision_cap=20):
        r"""
        Gives the kodaira symbol of the reduction at a given prime.

        INPUT:

        - ``prime`` -- A prime, i.e. a prime number if this Frey
          curve is defined over a subring of QQ or a maximal ideal
          of the ring of integers if it is defined over a subring
          of a number field.
        - ``verbose`` -- A boolean value or an integer
          (default: False). When set to True or any value
          larger then zero will print comments to stdout
          about the computations being done whilst busy. If
          set to False or 0 will not print such comments.
          If set to any negative value will also prevent
          the printing of any warnings.
          If this method calls any method that accepts an
          argument verbose will pass this argument to it.
          If such a method fulfills a minor task within
          this method and the argument verbose was larger
          than 0, will instead pass 1 less than the given
          argument. This makes it so a higher value will
          print more details about the computation than a
          lower one.
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
                           verbose=(verbose-1 if verbose>0 else verbose))
        result = performTatesAlgorithm(self,
                                       initial_values=Tp,
                                       coefficient_ring=self.base(),
                                       pAdics=pAdics,
                                       verbose=verbose,
                                       precision_cap=precision_cap,
                                       only_calculate=['minimal_model'])
        if len(result) == 1:
            return result[0][0][0]
        else:
            return ConditionalValue([(val[0], con) for val, con in result])
    
    @cached_method(key=lambda self, prime, verbose, precision_cap: (prime, precision_cap))
    def conductor_exponent(self, prime, verbose=False, precision_cap=20):
        r"""
        Gives the conductor exponent at a given prime.

        INPUT:

        - ``prime`` -- A prime, i.e. a prime number if this Frey
          curve is defined over a subring of QQ or a maximal ideal
          of the ring of integers if it is defined over a subring
          of a number field.
        - ``verbose`` -- A boolean value or an integer
          (default: False). When set to True or any value
          larger then zero will print comments to stdout
          about the computations being done whilst busy. If
          set to False or 0 will not print such comments.
          If set to any negative value will also prevent
          the printing of any warnings.
          If this method calls any method that accepts an
          argument verbose will pass this argument to it.
          If such a method fulfills a minor task within
          this method and the argument verbose was larger
          than 0, will instead pass 1 less than the given
          argument. This makes it so a higher value will
          print more details about the computation than a
          lower one.
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
            if len(result) == 1:
                return result[0]
            else:
                return ConditionalValue(result)
        
        pAdics = pAdicBase(self.definition_ring(), prime)
        Tp = self._initial_tree(pAdics.prime_below(self._R),
                                verbose=(verbose-1 if verbose>0 else verbose))
        result = performTatesAlgorithm(self,
                                       initial_values=Tp,
                                       coefficient_ring=self.base(),
                                       pAdics=pAdics,
                                       verbose=verbose,
                                       precision_cap=precision_cap,
                                       only_calculate=['conductor'])
        if(len(result) == 1):
            return result[0][0][0]
        else:
            return ConditionalValue([(val[0], con) for val,con in result])

    def base_extend(self, R):
        if (hasattr(R, 'domain') and
            R.domain() == self.definition_ring()):
            dom = self.base_ring()
            codom = dom.change_ring(R.codomain())
            F = RingHomomorphism_from_base(dom.Hom(codom), R)
            result = EllipticCurve_generic.base_extend(self, F)
        else:
            result = EllipticCurve_generic.base_extend(self, R)
        if ( (is_PolynomialRing(result.base_ring()) or
              is_MPolynomialRing(result.base_ring())) and
             (result.base_ring().variable_names() ==
              tuple(str(v) for v in self.parameters()))):
            return FreyCurve(result,
                             parameter_ring=self._R,
                             conversion=self._R_to_base,
                             condition=self._condition)
        return result

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
        - ``verbose`` -- A boolean value or an integer
          (default: False). When set to True or any value
          larger then zero will print comments to stdout
          about the computations being done whilst busy. If
          set to False or 0 will not print such comments.
          If set to any negative value will also prevent
          the printing of any warnings.
          If this method calls any method that accepts an
          argument verbose will pass this argument to it.
          If such a method fulfills a minor task within
          this method and the argument verbose was larger
          than 0, will instead pass 1 less than the given
          argument. This makes it so a higher value will
          print more details about the computation than a
          lower one.
        - ``precision_cap`` A strictly positive integer
          (default: 20) giving the maximal precision level to be
          used in p-Adic arithmetic.

        OUTPUT:

        A conditional expression that gives the conductor of
        this Frey curve for each possible value of the
        parameters on which it depends. The left side of this
        expression is some expression that gives the part of
        the conductor at all the primes given in additive_primes,
        whilst the right side is a string describing how to
        compute the part of the conductor that is coprime to
        those primes. The latter contains the operator Rad_P
        which refers to taking the radical of an expression
        ignoring those primes in additive_primes.
        """
        if additive_primes is None:
            additive_primes = self.primes_of_possible_additive_reduction()
        result = product(P^self.conductor_exponent(P, verbose=verbose,
                                                 precision_cap=precision_cap)
                         for P in additive_primes)
        return ConditionalExpression(['*','\cdot',2], result, "Rad_P( " + str(self.discriminant().factor()) + " )")
        
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

    def _repr_(self):
        """
        String representation of a Frey curve.

        REMARK:

        This is a direct copy from the code included
        in EllipticCurve_number_field
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
        s += "over %s "%(self.definition_ring(),)
        s += "with parameters %s"%(self.parameters(),)
        return s

class FreyQcurve(FreyCurve, Qcurve):
    r"""
    A Frey curve that is a Q-curve over some number field
    """
    def __init__(self, curve, parameter_ring=ZZ, conversion=None, condition=None, isogenies={}):
        r"""
        Initializes a Frey-Q-curve.

        This initialization calls the initialization of both
        the Qcurve and the FreyCurve class. Note however
        that for the Qcurve class the parameter guessed degrees
        is always set to zero as there is no good way to guess
        isogenies of a Frey curve of a given degree.

        The method _init_curve inside the class Qcurve is
        overwritten by this class's method, hence when the
        initialization of Qcurve is called, this method is
        called instead. This tricks the Qcurve class into
        thinking that this class is in fact defined over
        a number field even though it is not.

        INPUT:

        - ``curve`` -- An elliptic curve defined over some
          polynomial ring over a number field or any argument that
          would produce such a curve when passed to the
          constructor EllipticCurve. This curve will be taken
          over a polynomial ring over the minimal galois
          extension of its base field and will become the
          Frey Q-curve.
        - ``parameter_ring`` -- A ring (default: ZZ) that
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
         - ``isogenies`` -- A dictionary (default: {}) with as keys elements
           of the galois group of the base field of the Q-curve and as values
           data of the corresponding isogeny from the galois conjugate of this
           Q-curve to itself. This data can be either an isogeny as a Sage
           object or a tuple of an algebraic integer (defined as an element of
           some number field) and a strictly positive integer, which are
           respectively the $\lambda$ such that the isogeny is
           $z \mapsto \lambda z$ on the complex numbers and the degree of the
           isogeny.
        """
        FreyCurve.__init__(self, curve, parameter_ring=parameter_ring,
                           conversion=conversion, condition=condition)
        Qcurve.__init__(self, curve, isogenies=isogenies)

    def _init_curve(self, curve):
        r"""
        Initializes the curve of this Frey Q-curve.

        This overwrites the method _init_curve inside the
        class Qcurve. Note that most things have already
        been initialized by the class FreyCurve at this point,
        but we initialize again to make sure the following
        conditions are satisfied:
        - The curve is defined over a polynomial ring over some
          number field.
        - The number field is a galois field. If not it will be
          replaced by its galois closure.
        """
        K = self.definition_field()
        if not is_NumberField(K):
            raise ValueError("The ring %s is not a number field."%(K,))
        if not K.is_galois():
            Kgal = K.galois_closure(names=K.variable_name() + 'g')
            iota = K.hom([a.minpoly().change_ring(Kgal).roots()[0][0] for a in K.gens()], Kgal)
            ainvs = [a.change_ring(iota) for a in curve.a_invariants]
            S = self.base().change_ring(Kgal)
            EllipticCurve_generic.__init__(self, S, ainvs)

    def definition_field(self):
        r"""
        Gives the field over which this Frey Q-curve is defined.

        OUTPUT:

        The number field over which this Frey Q-curve is defined.
        """
        return self.definition_ring()

    def base_extend(self, R):
        result = FreyCurve.base_extend(self, R)
        if (isinstance(result, FreyCurve) and
            is_NumberField(result.definition_ring())):
            K = self.definition_field()
            L = result.definition_ring()
            r = K.gen().minpoly().change_ring(L).roots()
            if len(r) > 0:
                return FreyQcurve(result,
                                  isogenies=self._isogeny_data(L),
                                  parameter_ring=result._R,
                                  conversion=result._R_to_base,
                                  condition=result._condition)
        return result

    def twist(self, gamma):
        r"""
        Gives the twist of this Frey Q-curve by a given element gamma.

        INPUT:

        - ``gamma`` -- An element of a number field.

        OUTPUT:
        
        A Frey Q-curve defined over the composite field of the field over
        which this Frey Q-curve is completely defined and the parent of
        gamma, that is the twist of this Q-curve by gamma, i.e. if this
        Frey Q-curve was given by
        
        .. MATH::

        E : y^2 = x^3 + a_2 x^2 + a_4 x + a_6

        the twisted Q-curve is given by

        .. MATH::
        
        E : y^2 = x^3 + \gamma a_2 x^2 + \gamma^2 a_4 x + \gamma^3 a_6
        
        """
        K_E = self.complete_definition_field()
        K_gamma = gamma.parent()
        K, iota, gamma_map = composite_field(K_E, K_gamma, give_maps=True)
        gamma = gamma_map(gamma)
        E_map = iota * self._to_Kl
        E = twist_elliptic_curve(self.change_ring(E_map), gamma)
        ainvs = E.a_invariants()
        l = self.isogeny_lambda
        d = self.degree_map
        G = K.galois_group()
        isogenies = dict()
        for s in G:
            L, K_to_L, alpha = field_with_root(K, s(gamma)/gamma, give_embedding=True)
            isogenies[s] = (K_to_L(iota(l(s))) * alpha, d(s))
        H = [t for t in G if (all(t(isogenies[s][0]) == isogenies[s][0] for s in G) and
                              all(t(c) == c for a in ainvs for c in a.coefficients()))]
        Kmin = fixed_field(H)
        if Kmin != K:
            isogenies_min = {}
            G = Kmin.galois_group()
            for s in Kmin.galois_group():
                l, d = isogenies[galois_field_change(s, K)]
                isogenies_min[s] = (Kmin(l), d)
            ainvs = {a.change_ring(Kmin) for a in ainvs}
            conversion = (K_E.hom([a.minpoly().change_ring(Kmin).roots()[0][0]
                                  for a in K_E.gens()], Kmin) *
                          self_to_Kl * self._R_to_base)
            return FreyQcurve(ainvs, isogenies=isogenies_min,
                              parameter_ring=self._R,
                              conversion=conversion,
                              condition=self._condition)
        conversion = E_map * self._condition
        return FreyQcurve(ainvs, isogenies=isogenies,
                          parameter_ring=self._R,
                          conversion=self._R_to_base,
                          condition=self._condition)

    @cached_method
    def conductor_restriction_of_scalars(self, additive_primes=None,
                                         verbose=False, precision_cap=20):
        r"""
        Gives the conductor of the restriction of scalars of this Frey Q-curve.

        Note that since this is a Frey curve, the solution might depend on the
        parameter so the outcome will be a conditional expression that contains
        a factor expressed as a string since it can not be explicitly computed.

        INPUT:

        - ``bad_primes`` -- An iterable containing prime ideals
          or prime numbers, if the decomposition field is QQ, that
          contains all the primes at which this curve, over the
          decomposition field, can have additive reduction and all
          the primes at which the decomposition field ramifies. If
          set to None will compute this by using the method
          primes_of_possible_additive_reduction and by computing
          the ramified primes of the decomposition field.
        - ``verbose`` -- A boolean value or an integer
          (default: False). When set to True or any value
          larger then zero will print comments to stdout
          about the computations being done whilst busy. If
          set to False or 0 will not print such comments.
          If set to any negative value will also prevent
          the printing of any warnings.
          If this method calls any method that accepts an
          argument verbose will pass this argument to it.
          If such a method fulfills a minor task within
          this method and the argument verbose was larger
          than 0, will instead pass 1 less than the given
          argument. This makes it so a higher value will
          print more details about the computation than a
          lower one.
        - ``precision_cap`` A strictly positive integer
          (default: 20) giving the maximal precision level to be
          used in p-Adic arithmetic.

        OUTPUT:

        The conductor of the restriction of scalars of this curve over the
        decomposition field. This will be a conditional expression containing
        on the left side a (conditional) expression of the part of the
        conductor coming from primes in additive_primes, whilst the right hand
        side is a string describing how to compute the part of the conductor
        coming from primes coprime to the primes in additivie_primes. The
        latter contains the operator Rad_P which refers to taking the radical
        of an expression ignoring those primes in additive_primes.
        """
        K0 = self.definition_field()
        K = self.decomposition_field()
        if K0 != K:
            iota = K0.hom([a.minpoly().change_ring(K).roots()[0][0] for a in K0.gens()], K)
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
                        verbose=verbose,
                        precision_cap=precision_cap)
        additive_part = N.left()
        Dsqr = K.discriminant()^2
        if isinstance(additive_part, ConditionalExpression):
            additive_factors = N.left().factors()
            left_factors = {p: e * factors[f] for f in factors for p,e in f.absolute_norm().factor()}
            for p, e in Dsqr.factor():
                if p in left_factors:
                    left_factors[p] = left_factors[p] + e
                else:
                    left_factors[p] = e
            left = Dsqr.unit() * product(p^e for p,e in left_factors.iteritems())
        else:
            left = additive_part.absolute_norm() * Dsqr
        return ConditionalExpression(N.operator(),
                                     left,
                                     "Norm(" + N.right() +")")

    def _newform_levels(self, N=None, **kwds):
        r"""
        Gives the possible levels of newforms associated to this Frey Q-curve.

        NOTE:

        These newforms are computed using only finitely many information which can
        be supplied by the user through the parameter N. Otherwise the information
        will be limited to the left side of the expression returned by
        conductor_restriction_of scalars.

        For Frey curves this is sufficient information, since level lowering results
        often imply that a newform of the level mentioned above should exist. Note
        that this is however not checked nor proved by this code, hence the user
        should verify this themself and supply a different N if necessary.

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
        - ``N`` -- A non-negative integer (default: None) or ConditionalValue with
          such values giving the conductor of the restriction of scalars of this
          Q-curve over the decomposition field. Will be computed using the
          corresponding method if set to None, in which case only the non-radical
          part of the resulting ConditionalExpression will be used.
        
        OUTPUT:

        A list of tuples, each tuple representing one of the options for the levels
        of the newforms associated to this Q-curve. If a prime was given, these
        tuples will contain the respective exponent of the given prime for each
        newform. If no prime was given, they will contain the respective level of
        each newform.

        If the argument N was a ConditionalValue will return a list of tuples
        each of which contain a list as described above as their first entry
        and a condition for which this list is relevant as their second entry.
        """
        if N is None:
            N = self.conductor_restriction_of_scalars().left()
        if isinstance(N, ConditionalValue):
            return [(self._newform_levels(N=Ni, **kwds), con) for Ni, con in N]
        return Qcurve._newform_levels(self, N=N, **kwds)
        
    def _repr_(self):
        """
        String representation of a Frey Q-curve.

        REMARK:

        This is a direct copy from the code included
        in EllipticCurve_number_field
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
        s += "over %s "%(self.definition_ring(),)
        s += "with parameters %s"%(self.parameters(),)
        return s
