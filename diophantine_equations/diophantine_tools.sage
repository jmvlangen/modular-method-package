from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic

class Restriction(SageObject):
    """
    A class that represents a restriction on some variables.
    """
    
    def __init__(self, variables, exclude_primes=[]):
        self._vars = list(variables)
        if hasattr(exclude_primes, '__iter__'):
            self._excluded = []
            for p in exclude_primes:
                if p in ZZ.ideal_monoid():
                    p = p.gens()[0]
                self._excluded.append(p)
        else:
            raise ValueError("The primes to be excluded are not in a list.")
        
    def variables(self):
        return self._vars
        
    def _repr_(self):
        return "A restriction on the variables %s."%(self._vars,)
        
    def should_exclude(self, pAdics):
        K = pAdics.number_field()
        P = pAdics.prime_ideal()
        for Q in self._excluded:
            if Q in ZZ or Q.number_field().is_subring(K):
                if K is QQ:
                    if P.gens()[0] == Q:
                        return True
                else:
                    if P in K.primes_above(Q):
                        return True
            if K.is_subring(Q.number_field()):
                if Q in Q.number_field().primes_above(P):
                    return True
        return False
        
    def apply_to_tree(self, pAdic_tree, **kwds):
        """
        Applies this restriction to a pAdic_tree.
        
        Additional keywords can be defined by individual instances
        of this function. Therefore all remaining keywords have
        to be captured using `**kwds`.
        """
        raise NotImplementedError("The method 'apply_to_tree' was not (correctly) implemented for this restriction.")
        
class PolynomialRestriction(Restriction):
    """
    A polynomial restriction on some variables
    """
    
    def __init__(self, polynomial, eq_zero=True, exclude_primes=[]):
        if not isinstance(polynomial, Polynomial) and \
           not isinstance(polynomial, MPolynomial):
            raise ValueError("The given argument %s is not a polynomial."%(polynomial,))
        self._f = polynomial
        self._eq_zero = eq_zero
        Restriction.__init__(self, self._f.variables(),
                             exclude_primes=exclude_primes)
        
    def polynomial(self):
        return self._f
        
    def _repr_(self):
        if self._eq_zero:
            return "The restriction %s == 0 on the variables %s."%(self.polynomial(),
                                                                   self.variables())
        else:
            return "The restriction %s != 0 on the variables %s."%(self.polynomial(),
                                                                   self.variables())
                                                                   
    def _latex_(self):
        return latex(polynomial) + ("=" if self._eq_zero else "\\neq") + \
               latex(0)
               
    def apply_to_tree(self, pAdic_tree, pAdics=None, precision=20,
                      verbose=False, **kwds):
        if self.should_exclude(pAdics):
            return
        if "precision_cap" in kwds:
            precision = kwds[precision_cap]
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        if self._eq_zero:
            pAdic_tree.apply_polynomial_restriction(self._f,
                                                    pAdics=pAdics,
                                                    precision=precision,
                                                    verbose=verbose)
        else:
            pAdic_tree.apply_polynomial_anti_restriction(self._f,
                                                         pAdics=pAdics,
                                                         precision=precision,
                                                         verbose=verbose)
        
class ComparisonRestriction(Restriction):
    """
    A comparison restriction on some variables
    
    For one variable this will act as the usual comparison on that variable.
    For multiple variables, the following rule will be applied:
        (a_1, ... , a_n) ~ (b_1, ... , b_n) if and
        only if a_i ~ b_i for all i,
    where ~ can stand for any relation supported by this class.
    """
    
    EQUAL = 0
    GREATER = 2
    GREATER_EQUAL = 1
    LESSER = -2
    LESSER_EQUAL = -1
    
    def __init__(self, variables, comparison, values, exclude_primes=[]):
        if len(variables) != len(values):
            raise ValueError("The amount of variables must correspond to the amount of values.")
        if comparison not in [EQUAL, GREATER, GREATER_EQUAL,
                              LESSER, LESSER_EQUAL]:
            raise ValueError("%s is not a valid comparison."%(comparision,))
        Restriction.__init__(self, variables, exclude_primes=exclude_primes)
        self._r = comparison
        self._vals = values
        
    def _repr_(self):
        if self._r == EQUAL:
            return "The restriction %s == %s"%(tuple(self.variables()),
                                               tuple(self._vals))
        elif self._r == GREATER:
            return "The restriction %s > %s"%(tuple(self.variables()),
                                              tuple(self._vals))
        elif self._r == GREATER_EQUAL:
            return "The restriction %s >= %s"%(tuple(self.variables()),
                                               tuple(self._vals))
        elif self._r == LESSER:
            return "The restriction %s < %s"%(tuple(self.variables()),
                                              tuple(self._vals))
        elif self._r == LESSER_EQUAL:
            return "The restriction %s <= %s"%(tuple(self.variables()),
                                               tuple(self._vals))
        else:
            raise ValueError("%s is not a valid comparison."%(self._r,))
            
    def _latex_(self):
        if self._r == EQUAL:
            return latex(tuple(self.variables())) + \
                   "=" + latex(tuple(self._vals))
        elif self._r == GREATER:
            return latex(tuple(self.variables())) + \
                   ">" + latex(tuple(self._vals))
        elif self._r == GREATER_EQUAL:
            return latex(tuple(self.variables())) + \
                   "\\geq" + latex(tuple(self._vals))
        elif self._r == LESSER:
            return latex(tuple(self.variables())) + \
                   "<" + latex(tuple(self._vals))
        elif self._r == LESSER_EQUAL:
            return latex(tuple(self.variables())) + \
                   "\\leq" + latex(tuple(self._vals))
        else:
            raise ValueError("%s is not a valid comparison."%(self._r,))
            
class CongruenceRestriction(Restriction):
    """
    A congruence restriction on some variables.
    """
    
    def __init__(self, variable, values, modulus, eq=True, exclude_primes=[]):
        if not hasattr(values, '__iter__'):
            values = [values]
        values = list(values)
        Restriction.__init__(self, [variable], exclude_primes=exclude_primes)
        self._vals = values
        self._mod = modulus
        self._eq = eq
        
    def _repr_(self):
        if self._eq:
            return "The restriction %s == %s modulo %s."%(self.variables()[0],
                                                          self._vals,
                                                          self._mod)
        else:
            return "The restriction %s != %s modulo %s."%(self.variables()[0],
                                                          self._vals,
                                                          self._mod)
    def _latex_vars(variables, values, symbol):
        result = ""
        n = len(variables)
        for i in range(n):
            result += latex(variables[i])
            result += " " + symbol + " "
            result += latex(values[i])
            if i < n - 1:
                result += " \\wedge "
        return result
                           
    def _latex_(self):
        result = ""
        symbol = ("\\equiv" if self._eq else "\\not\\equiv")
        n = len(values)
        if n > 1:
            result += "\\left\\{ \\begin{array}{ll}"
        for i in range(n):
            result += _latex_vars(variables, values[i], symbol)
            result += " & "
            result += "\\text{(mod }" + latex(self._mod) + ")"
            if i < n - 1:
                result += " \\\n"
        if n > 1:
            result += " \\end{array} \\right."
        return result
        
    def apply_to_tree(self, pAdic_tree, **kwds):
        if self.should_exclude(pAdic_tree.pAdics()):
            return
        tree_vars = pAdic_tree.variables()
        var = self._vars[0]
        if var in tree_vars:
            if self._eq:
                pAdic_tree.apply_congruence_restriction(var,
                                                        self._vals,
                                                        self._mod)
            else:
                pAdic_tree.apply_congruence_anti_restriction(var,
                                                             self._vals,
                                                             self._mod)
                                                            
class PolynomialCongruenceRestriction(Restriction):
    """
    A congruence restriction on some polynomial.
    """
    
    def __init__(self, polynomial, modulus, eq=True, exclude_primes=[]):
        if not isinstance(polynomial, Polynomial) and \
           not isinstance(polynomial, MPolynomial):
            raise ValueError("The given argument %s is not a polynomial."%(polynomial,))
        self._f = polynomial
        Restriction.__init__(self, self._f.variables(),
                             exclude_primes=exclude_primes)
        self._mod = modulus
        self._eq = eq
        
    def _repr_(self):
        if self._eq:
            return "The restriction %s == 0 modulo %s."%(self._f,
                                                          self._mod)
        else:
            return "The restriction %s != 0 modulo %s."%(self._f,
                                                          self._mod)
    def _latex_vars(variables, values, symbol):
        result = ""
        n = len(variables)
        for i in range(n):
            result += latex(variables[i])
            result += " " + symbol + " "
            result += latex(values[i])
            if i < n - 1:
                result += " \\wedge "
        return result
                           
    def _latex_(self):
        result = ""
        symbol = ("\\equiv" if self._eq else "\\not\\equiv")
        n = len(values)
        if n > 1:
            result += "\\left\\{ \\begin{array}{ll}"
        for i in range(n):
            result += _latex_vars(variables, values[i], symbol)
            result += " & "
            result += "\\text{(mod }" + latex(self._mod) + ")"
            if i < n - 1:
                result += " \\\n"
        if n > 1:
            result += " \\end{array} \\right."
        return result
        
    def apply_to_tree(self, pAdic_tree, pAdics=None, verbose=False, **kwds):
        if self.should_exclude(pAdics):
            return
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        n = pAdics.valuation(self._mod)
        if n > 0:
            if self._eq:
                pAdic_tree.apply_polynomial_restriction(self._f,
                                                        pAdics=pAdics,
                                                        precision=n,
                                                        verbose=verbose)
            else:
                pAdic_tree.apply_polynomial_anti_restriction(self._f,
                                                             pAdics=pAdics,
                                                             precision=n,
                                                             verbose=verbose)

class CoprimeRestriction(Restriction):

    def __init__(self, variables, coprimality=2, exclude_primes=[]):
        Restriction.__init__(self, variables, exclude_primes=exclude_primes)
        if coprimality not in ZZ or \
           coprimality < 2 or \
           coprimality > len(self._vars):
            raise ValueError("The coprimality should be an integer between 2 and the number of variables, not %s."%(coprimality,))
        self._coprime = coprimality
        
    def _repr_(self):
        return "Restriction making the variables " + \
               str(self._vars) + " " + \
               (str(self._coprime) + "-wise" if self._coprime > 2 else "pairwise") + \
               " coprime."
               
    def apply_to_tree(self, pAdicTree, **kwds):
        if self.should_exclude(pAdicTree.pAdics()):
            return
        pAdicTree.apply_coprimality_restriction(variables=self._vars,
                                                coPrimality=self._coprime)
                                                                
class PolynomialPowerRestriction(Restriction):
    """
    A restriction on variables of the form:
        f = z^n
    where f is a polynomial in some variables except z,
    z is another (unknown) variable,
    and n is a parameter with at least a certain strictly positive value.
    """
    
    def __init__(self, polynomial, least_exponent, power_variable=None,
                 exponent_variable=None, exclude_primes=[]):
        if not isinstance(polynomial, Polynomial) and \
           not isinstance(polynomial, MPolynomial):
            raise ValueError("The given argument %s is not a polynomial."%(polynomial,))
        self._f = polynomial
        self._pow_var = power_variable
        if least_exponent not in ZZ or least_exponent <= 0:
            raise ValueError("The least exponent should be a strictly positive integer, not %s."%(least_exponent,))
        self._exp = least_exponent
        self._exp_var = exponent_variable
        variables = self._f.variables()
        if not self._pow_var is None:
            variables.append(self._pow_var)
        Restriction.__init__(self, variables, exclude_primes=exclude_primes)
        
    def _repr_(self):
        return str(self._f) + \
               " = " + \
               (str(self._pow_var) if not self._pow_var is None else "z") + \
               "^(" + \
               (str(self._exp_var) if not self._exp_var is None else "n") + \
               ")"
               
    def _latex_(self):
        return latex(self._f) + \
               " = " + \
               (latex(self._pow_var) if not self._pow_var is None else "z") + \
               "^{" + \
               (latex(self._exp_var) if not self._exp_var is None else "n") + \
               "}"
               
    def apply_to_tree(self, pAdicTree, pAdics=None, verbose=False, **kwds):
        if self.should_exclude(pAdics):
            return
        orig_vars = pAdicTree.variables()
        variables = list(Set(orig_vars).union(Set(self.variables())))
        pAdicTree.change_variables_to(variables)
        if pAdics is None or self._pow_var in variables:
            pAdics = pAdicTree.pAdics()
        R = PolynomialRing(pAdics.number_field(), variables);
        f = R(self._f)
        T1, T2 = find_pAdicRoots(f,
                                 pAdics=pAdics,
                                 variables=variables,
                                 value_tree=pAdicTree,
                                 precision=1,
                                 verbose=verbose)
        T3, T4 = find_pAdicRoots(f,
                                 pAdics=pAdics,
                                 variables=variables,
                                 value_tree=pAdicTree,
                                 precision=self._exp,
                                 verbose=verbose)
        T2.merge(T3)
        pAdicTree._root = T2
        pAdicTree.change_variables_to(orig_vars)
   
def get_variables_in_polynomials(polynomial_list):
    variables = Set([])
    for poly in polynomial_list:
        variables = variables.union(Set(poly.variables()))
    return list(variables)

def _short_ideal_string(ideal):
    if ideal.is_principal():
        s = str(ideal.gens_reduced()[0])
        if s.count('+') > 0 or s.count('-') > 0 \
                            or s.count('*') > 0 \
                            or s.count('^') > 0:
            s = '(' + s + ')'
        return s
    else:
        return ideal._repr_short()
        
def _short_tuple_string(t):
    if len(t) == 1:
        return str(t[0])
    else:
        return str(t)

def _gcd_coefficients(polynomial):
    R = polynomial.base_ring()
    if R is QQ:
        R = NumberField(PolynomialRing(QQ,'x').gen())
    I = sum([R.ideal(c) for c in polynomial.coefficients()])
    return I
  
class ConditionalConductor(SageObject):
    
    def __init__(self, conditions, radical_polynomial, variables):
        self._dict = conditions
        self._rad = radical_polynomial
        self._vars = tuple(variables)
        self._reduce_coefficient()
        self._P = []
        self._conditions = []
        self._constants = []       
        ls = list(self._dict.iteritems())
        self._str_init(ls)
    
    def _reduce_coefficient(self):
        I = _gcd_coefficients(self._rad)
        if I == 0:
            return
        if self._rad.base_ring().degree() == 1:
            I = QQ(I.gens_reduced()[0])
        J0 = I^0
        for (p,n) in I.factor():
            if n > 0:
                J0 *= p^n
                if p in ZZ:
                    p = ZZ.ideal(p)
                if p not in self._dict:
                    self._dict[p] = [[1, pAdicTree(variables=self._vars,
                                                   prime=p)]]
        J = self._rad.base_ring().ideal(J0)
        while not J.is_principal():
            J *= J0
        k = J.gens_reduced()[0]
        self._rad *= (1/k)
    
    def _init_constant(self, p, con, e_index):
        if len(con) <= 1 and con[0][0] == 0:
            return
        result = _short_ideal_string(p)
        if len(con)> 1 or con[0][0] != 1:
            result += "^("
            result += ("e" + str(e_index) if len(con) > 1 else str(con[0][0]))
            result += ")"
        self._constants.append(result)
        
    def _init_condition(self, p, con, e_index):
        if len(con) > 1:
            result = ""
            label = "with e" + str(e_index) + " = "
            result += label
            label_white = "\n" + (" " * len(label))
            for j in range(len(con)):
                c = con[j]
                if not c[1].root().is_full():
                    vals, mod = c[1].give_as_congruence_condition()
                    if j > 0:
                        result += label_white
                    result += str(c[0]) + " if "
                    result += _short_tuple_string(self._vars) + " = "
                    for i in range(len(vals)):
                        result += _short_tuple_string(vals[i])
                        if i < len(vals) - 1:
                            result += ", "
                    result += " (mod " + _short_ideal_string(mod) + ")"
            self._conditions.append(result)
    
    def _str_init(self, ls):
        e_index = 1
        for (p, con) in ls:
            self._P.append(p)
            self._init_constant(p, con, e_index)
            self._init_condition(p, con, e_index)
            if len(con) > 1:
                e_index += 1
        
    def _conductor_string(self):
        result = ""
        for c in self._constants:
            result += c
            result += " * "
        result += "Rad_P(" + str(self._rad.factor()) + ")\n"
        result += "where P = {"
        for i in range(len(self._P)):
            result += _short_ideal_string(self._P[i])
            if i < len(self._P) - 1:
                result += ", "
        result += "}\n"
        return result
        
    def _condition_string(self):
        result = ""
        for c in self._conditions:
            result += c + "\n"
        return result
                
    def _repr_(self):
        return self._conductor_string() + self._condition_string()
    
class DiophantineAnalyzer(SageObject):
    
    def __init__(self, variables, parameters=[], value_ring=ZZ):
        if not hasattr(variables, "__iter__"):
            variables = [variables]
        self._poly_ring = PolynomialRing(value_ring, variables)
        field_poly_ring = PolynomialRing(self.variable_field(),
                                         self.variables())
        field_poly_ring.inject_variables()
        if not hasattr(parameters, "__iter__"):
            parameters = [parameters]
        self._par = parameters
        self._restr = []
        self._tree_dict = {}
        
    def variables(self):
        return self._poly_ring.gens()
        
    def variable_ring(self):
        return self._poly_ring.base()
        
    def variable_field(self):
        return self.variable_ring().fraction_field()
        
    def parameters(self):
        return tuple(self._par)
        
    def add_restriction(self, restriction):
        if not isinstance(restriction, Restriction):
            raise ValueError("The given restriction %s is not a restriction."%(restriction,))
        self._restr.append(restriction)
        self._tree_dict = {}
        
    def _is_prime(self, prime):
        return prime in self.variable_field().ideal_monoid() and prime.is_prime()
        
    def get_initial_tree(self, pAdics, variables=None, verbose=False,
                         precision_cap=20, **kwds):
        prime = pAdics.prime_ideal()
        if pAdics.number_field() is QQ:
            prime = prime.gens()[0]
        key = (prime if variables is None else (prime, tuple(variables)))
        if key not in self._tree_dict:
            if variables is None or prime not in self._tree_dict:
                T = pAdicTree(variables=self.variables(),
                              prime=self._get_prime_below(prime,
                                                          pAdics.number_field()))
                for restriction in self._restr:
                    restriction.apply_to_tree(T, pAdics=pAdics, 
                                              precision_cap=precision_cap,
                                              verbose=verbose, **kwds)
                self._tree_dict[prime] = T
            if not variables is None:
                T = self._tree_dict[prime].copy()
                T.change_variables_to(variables)
                self._tree_dict[key] = T
        return self._tree_dict[key].copy()
    
    def _determine_special_primes(self, K, c4, Delta, variables):
        if len(variables) == 1:
            result = c4.resultant(Delta)
            return K.ideal(result).prime_factors()
        elif len(variables) == 2 and c4.is_homogeneous() and Delta.is_homogeneous():
            result = c4.macaulay_resultant(Delta)
            return K.ideal(result).prime_factors()
        else:
            I = sum([K.ideal(c) for c in c4.coefficients()])
            J = sum([K.ideal(c) for c in Delta.coefficients()])
            return (I + J).prime_factors()
            
    def _get_prime_below(self, P, L):
        K = self.variable_field()
        if not K.is_subring(L):
            raise ValueError("%s is not an extension of %s."%(L, K))
        if L is QQ:
            return P
        else:
            p = P.smallest_integer()
            if K is QQ:
                return p
            else:
                for Q in K.primes_above(p):
                    if P in L.primes_above(Q):
                        return P
                raise ValueError("No prime lies below %s."%(P))
      
    def compute_conductor(self, E, precision_cap=20, verbose=False, model=False):
        if not isinstance(E, EllipticCurve_generic):
            raise ValueError("%s is not an elliptic curve."%(E,))
        variables = get_variables_in_polynomials(E.a_invariants())
        for var in variables:
            if var not in self._poly_ring:
                raise ValueError("The variable %s is not part of this analyzer."%(var,))
        K = E.base().base()
        if K is QQ:
            R = PolynomialRing(QQ,'x')
            K.<a> = NumberField(R.gen())
        ring = PolynomialRing(K, variables)
        c4 = ring(E.c4())
        Delta = ring(E.discriminant())
        P = self._determine_special_primes(K, c4, Delta, variables)
        K = E.base().base()
        ring = PolynomialRing(K, variables)
        if K is QQ:
            P = [ZZ.ideal(p.gens()[0]) for p in P]
        print "Assuming \n\t%s\nand\n\t%s\nto be coprime outside\n%s"%(c4.factor(),
                                                                       Delta.factor(),
                                                                       P)
        conductor = {}
        minimal_model = {}
        for p in P:
            pAdics = pAdicBase(K, p)
            Tp = self.get_initial_tree(pAdics, variables,
                                       precision_cap=precision_cap,
                                       verbose=verbose)
            answer = performTatesAlgorithm(E, initial_values=Tp,
                                           coefficient_ring=ring,
                                           pAdics=pAdics,
                                           verbose=verbose,
                                           precision_cap=precision_cap,
                                           only_calculate=['conductor',
                                                           'minimal_model'])
            conductor[p] = [[x[0],x[2]] for x in answer]
            minimal_model[p] = [[x[1],x[2]] for x in answer]
        conductor_object = ConditionalConductor(conductor, Delta, variables)
        if model:
            return conductor_object, minimal_model
        else:
            return conductor_object
