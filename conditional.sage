class Condition(SageObject):
    """
    A class that represents a condition that must hold on some variables.
    """
    
    def __init__(self, variables):
       self._vars = list(variables)
        
    def variables(self):
        return self._vars
        
    def _repr_short(self):
        raise NotImplementedError("The method '_repr_short' was not (correctly) implemented for this condtion.")
        
    def _repr_(self):
        return "A condition on the variables %s."%(self._vars,)
        
    def apply_to_tree(self, pAdic_tree, **kwds):
        """
        Applies this condition to a pAdic_tree.
        
        Additional keywords can be defined by individual instances
        of this function. Therefore all remaining keywords have
        to be captured using `**kwds`.
        """
        raise NotImplementedError("The method 'apply_to_tree' was not (correctly) implemented for this condition.")
        
class PolynomialCondition(Condition):
    """
    A polynomial condition on some variables
    """
    
    def __init__(self, polynomial, eq_zero=True):
        if not isinstance(polynomial, Polynomial) and \
           not isinstance(polynomial, MPolynomial):
            raise ValueError("The given argument %s is not a polynomial."%(polynomial,))
        self._f = polynomial
        self._eq_zero = eq_zero
        Condition.__init__(self, self._f.gens())
        
    def polynomial(self):
        return self._f
        
    def _repr_short(self):
        return str(self.polynomial()) + \
               (" == " if self._eq_zero else " != ") + "0"
        
        
    def _repr_(self):
        return "The condition %s on the variables %s."%(self._repr_short,
                                                        self.variables())
                                                                   
    def _latex_(self):
        return latex(polynomial) + ("=" if self._eq_zero else "\\neq") + \
               latex(0)
               
    def apply_to_tree(self, pAdic_tree, precision=20, verbose=False, **kwds):
        if self._eq_zero:
            pAdic_tree.apply_polynomial_condition(self._f,
                                                    precision=precision,
                                                    verbose=verbose)
        else:
            pAdic_tree.apply_polynomial_anti_condition(self._f,
                                                         precision=precision,
                                                         verbose=verbose)
        
class ComparisonCondition(Condition):
    """
    A comparison condition on some variables
    
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
    
    def __init__(self, variables, comparison, values):
        if len(variables) != len(values):
            raise ValueError("The amount of variables must correspond to the amount of values.")
        if comparison not in [EQUAL, GREATER, GREATER_EQUAL,
                              LESSER, LESSER_EQUAL]:
            raise ValueError("%s is not a valid comparison."%(comparision,))
        Condition.__init__(self, variables)
        self._r = comparison
        self._vals = values
    
    def _repr_short(self):
        if self._r == EQUAL:
            return "%s == %s"%(tuple(self.variables()),
                               tuple(self._vals))
        elif self._r == GREATER:
            return "%s > %s"%(tuple(self.variables()),
                              tuple(self._vals))
        elif self._r == GREATER_EQUAL:
            return "%s >= %s"%(tuple(self.variables()),
                               tuple(self._vals))
        elif self._r == LESSER:
            return "%s < %s"%(tuple(self.variables()),
                              tuple(self._vals))
        elif self._r == LESSER_EQUAL:
            return "%s <= %s"%(tuple(self.variables()),
                               tuple(self._vals))
        else:
            raise ValueError("%s is not a valid comparison."%(self._r,))   
        
    def _repr_(self):
        return "The condition %s"%(self._repr_short()),
            
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
            
class CongruenceCondition(Condition):
    """
    A congruence condition on some variables.
    """
    
    def __init__(self, variable, values, modulus, eq=True):
        if not hasattr(values, '__iter__'):
            values = [values]
        values = list(values)
        Condition.__init__(self, [variable])
        self._vals = values
        self._mod = modulus
        self._eq = eq
        
    def _repr_short(self):
        result = str(self._vars)
        result += (" == " if self._eq else " != ")
        for i in range(len(self._vals)):
            result += str(self._vals[i])
            if i < len(self._vals) - 1:
                result += ", "
        result += " (mod "
        result += (self._mod._repr_short() if hasattr(self._mod, '_repr_short') else str(self._mod))
        result += ")"
        return result
        
    def _repr_(self):
        return "The condition %s on the variables %s."%(self._repr_short)
    
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
        
    def apply_to_tree(self, pAdicTree, **kwds):
        tree_vars = pAdicTree.variables()
        var = self._vars[0]
        if var in tree_vars:
            if self._eq:
                pAdicTree.apply_congruence_condition(var,
                                                       self._vals,
                                                       self._mod)
            else:
                pAdicTree.apply_congruence_anti_condition(var,
                                                            self._vals,
                                                            self._mod)

class CoprimeCondition(Condition):

    def __init__(self, variables, coprimality=2):
        Condition.__init__(self, variables)
        if coprimality not in ZZ or \
           coprimality < 2 or \
           coprimality > len(self._vars):
            raise ValueError("The coprimality should be an integer between 2 and the number of variables, not %s."%(coprimality,))
        self._coprime = coprimality
        
    def _repr_(self):
        return "Condition making the variables " + \
               str(self._vars) + " " + \
               (str(self._coprime) + "-wise" if self._coprime > 2 else "pairwise") + \
               " coprime."
               
    def apply_to_tree(self, pAdicTree, **kwds):
        pAdicTree.apply_coprimality_condition(variables=self._vars,
                                                coPrimality=self._coprime)
                                                                
class PolynomialPowerCondition(Condition):
    """
    A condition on variables of the form:
        f = z^n
    where f is a polynomial in some variables except z,
    z is another (unknown) variable,
    and n is a parameter with at least a certain strictly positive value.
    """
    
    def __init__(self, polynomial, least_exponent, power_variable=None,
                 exponent_variable=None):
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
        Condition.__init__(self, variables)
        
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
               
    def apply_to_tree(self, pAdicTree, verbose=False, **kwds):
        T1, T2 = find_pAdicRoots(self._f,
                                 pAdicTree.pAdics().prime_ideal(),
                                 variables=self.variables(),
                                 value_tree=pAdicTree,
                                 precision=1,
                                 verbose=verbose)
        T3, T4 = find_pAdicRoots(self._f,
                                 pAdicTree.pAdics().prime_ideal(),
                                 variables=self.variables(),
                                 value_tree=pAdicTree,
                                 precision=self._exp,
                                 verbose=verbose)
        T2.merge(T3)
        pAdicTree._root = T2
    
class ConditionalObject(SageObject):
    
    def __init__(self):
        raise NotImplementedError()
        
    def getValue(self):
        raise NotImplementedError()
    
    def get_value_condition_list(self):
        raise NotImplementedError()
        
    def _repr_(self):
        raise NotImplementedError()
        
    def _latex_(self):
        raise NotImplementedError()

class ConditionalValue(ConditionalObject):
    
    def __init__(self, value_condition_list):
        value_condition_list = list(value_condition_list)
        for (value, condition) in value_condition_list:
            if not isinstance(condition, Condition):
                raise ValueError("%s is not a valid condition."%(condition,))
        self._list = value_condition_list
    
    def get_value_condition_list(self):
        return self._list
    
    def _repr_(self):
        result = ""
        for i in range(len(self._list)):
            result += str(self._list[i][0])
            result += "\tif "
            result += self._list[i][1]._repr_short()
            if i < len(self._list) - 1:
                result += "\n"
        return result

    def _latex_(self):
        result = ""
        n = len(self._vals)
        if n > 1:
            result += "\\left\\{ \\begin{array}{ll}"
        for i in range(n):
            result += latex(self._vals[i][0])
            result += " & \\text{ if }"
            result += latex(self._vals[i][1])
        if n > 1:
            result += "\\end{array} \\right\\}"
        return result

def _value_condition_list(x):
    if instanceof(x, ConditionalObject):
        return x.get_value_condition_list
    else:
        return [(x, None)]
        
class Conditional_Power(ConditionalObject):
    
    def __init__(self, base, exponent):
        self._base = base
        self._exp = exponent
        
    def get_value_condition_list(self):
        result = []
        
