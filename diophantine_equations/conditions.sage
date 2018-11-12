from sage.rings.polynomial.multi_polynomial import is_MPolynomial
from sage.rings.polynomial.polynomial_element import is_Polynomial
import re
import itertools

class Condition_base(SageObject):
    """
    A class that represents a condition on some variables.
    """

    def __init__(self, variables):
        r"""
        The default constructor of a Condition.

        INPUT:
        
        - ``variables`` -- A collection of variables on
          which this condition applies. This may be any
          form of a variable, but will be converted into
          strings. Multiple variables with the same name
          are therefore not very well supported and may
          cause unpredictable behavior.
        """
        self._vars = tuple(str(v) for v in variables)

    def variables(self):
        r"""
        Gives the variables on which this Condition applies

        OUTPUT:

        A tuple of variables on which this Condition applies.
        """
        return self._vars

    def __and__(self, other):
        r"""
        Creates the condition that both conditions hold.

        INPUT:
        
        - ``other`` -- A Condition.

        OUTPUT:

        A Condition object that holds on all values where
        both this Condition object and the given Condition
        object hold.
        """
        return AndCondition(self, other)

    def __or__(self, other):
        r"""
        Creates the condition that either condition holds.

        INPUT:
        
        - ``other`` -- A Condition.

        OUTPUT:

        A Condition object that holds on all values where
        either this Condition object or the given Condition
        object holds.
        """
        return OrCondition(self, other)

    def __invert__(self):
        r"""
        Creates the condition that this condition does not hold.

        OUTPUT:

        A Condition object that holds on all values where
        this Condition does not hold.
        """
        return NotCondition(self)

    def _repr_(self):
        return "A condition on the variables %s"%(self._vars,)

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree that satisfies
        the condition defined by this object on the variables
        therein.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the
        given pAdicTree which satisfies the condition
        defined by this object. If complement was set to
        True will return a tuple with the afore mentioned
        as its first entry and the complement of that tree
        within the given pAdicTree as its second argument.
        """
        raise NotImplementedError("The method 'apply_to_tree' is not implemented for the base condition class")

class PolynomialCondition(Condition_base):
    r"""
    A condition that a certain polynomial is zero.
    """

    def __init__(self, polynomial):
        r"""
        The constructor of a PolynomialCondition.

        INPUT:

        - ``polynomial`` -- A polynomial of one or more variables
          which should be zero or non-zero
        """
        if not (is_Polynomial(polynomial) or is_MPolynomial(polynomial)):
            raise ValueError("The given argument %s is not a polynomial."%(polynomial,))
        self._f = polynomial
        Condition_base.__init__(self, self._f.variables())

    def polynomial(self):
        r"""
        OUTPUT:

        The polynomial that defines this PolynomialCondition.
        """
        return self._f

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   precision=20, verbose=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those
        values for the variable such that the polynomial
        of this condition is zero on them.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.
        - ``precision`` -- A strictly positive integer
          (default: 20) giving up to what precision the
          resulting tree should be found.
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

        A pAdicTree object that contatins that part of the
        given pAdicTree which satisfies the polynomial of
        this condition being equal to zero modulo the
        prime defined by pAdics up to the power precision.
        If complement is set to True will also give the
        complement of this tree in the given tree as a 
        second return value.
        """
        if "precision_cap" in kwds and precision > kwds["precision_cap"]:
            print "Warning: A p-Adic tree for %s is computed with a lower precision than the given precision %s, due to a given precision_cap %s!"%(self, precision, kwds["precision_cap"])
            precision = kwds["precision_cap"]
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics,
                                   full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        big_vars = list(pAdic_tree.variables())
        for var in self.variables():
            if var not in big_vars:
                big_vars.append(var)
        K = pAdics.number_field()
        K0 = self.polynomial().parent().base()
        if K0.is_subring(QQ):
            iota = K0.hom(K)
        else:
            iota = K0.hom([a.minpoly().change_ring(K).roots()[0][0] for a in K0.gens()], K)
        S = PolynomialRing(K, big_vars)
        Tyes, Tno = find_pAdicRoots(S(self.polynomial().change_ring(iota)),
                                    pAdics = pAdics,
                                    variables=[S(v) for v in pAdic_tree.variables()],
                                    value_tree=pAdic_tree.root(),
                                    precision=precision,
                                    verbose=verbose)
        Tyes = pAdicTree(variables=big_vars, root=Tyes)
        Tyes = Tyes.change_variables_to(pAdic_tree.variables())
        if complement:
            Tno = pAdicTree(variables=big_vars, root=Tno)
            Tno = Tno.change_variables_to(pAdic_tree.variables())
            return Tyes, Tno
        else:
            return Tyes

    def _repr_(self):
        return "The condition that %s == 0"%(self.polynomial(),)

    def _repr_short(self):
        return "%s == 0"%(self.polynomial(),)
        
    def _latex_(self):
        return latex(self.polynomial()) + " = " + latex(0)

    def _cache_key(self):
        return 'TreeCondition', self.polynomial()

class CongruenceCondition(PolynomialCondition):
    r"""
    Defines a congruence that should hold for a given polynomial.
    """

    def __init__(self, polynomial, modulus):
        r"""
        A constructor of a CongruenceCondition.

        INPUT:

        - ``polynomial`` -- A polynomial for which a
          congruence should hold.
        - ``modulus`` -- An algebraic integer or integral
          ideal of the ring of integers of a number field
          modulo which the polynomial should be zero.
        """
        PolynomialCondition.__init__(self, polynomial)
        self._mod = modulus
        if hasattr(self._mod, 'is_principal') and self._mod.is_principal():
            self._mod = self._mod.gens_reduced()[0]    

    def modulus(self):
        r"""
        OUTPUT:

        The algebraic integer or integral ideal modulo which
        the polynomial of this condition should be zero.
        """
        return self._mod
        
    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   verbose=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those
        values for the variables such that the polynomial
        of this condition is zero on them modulo the
        modulus given in this condition.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.
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

        A pAdicTree object that contatins that part of the
        given pAdicTree which satisfies the polynomial of
        this condition being equal to zero modulo the
        maximal prime power of the prime in the given
        pAdics that divides the modulus given in this
        Condition.
        If complement is set to True will also give the
        complement of this tree or the same tree if the
        prime in the given pAdics does not divide the
        modulus.
        """
        precision = pAdics.valuation(self.modulus())
        result = PolynomialCondition.pAdic_tree(self, pAdic_tree=pAdic_tree,
                                                pAdics=pAdics, complement=complement,
                                                verbose=verbose, precision=precision,
                                                **kwds)
        if complement and precision == 0:
            return result[0], result[0]
        else:
            return result
    
    def _repr_(self):
        mod = self.modulus()
        mod_str = (mod._repr_short() if hasattr(mod, '_repr_short') else str(mod))
        return "The condition that %s == 0 modulo %s"%(self.polynomial(),
                                                       mod_str)

    def _repr_short(self):
        mod = self.modulus()
        mod_str = (mod._repr_short() if hasattr(mod, '_repr_short') else str(mod))
        return "%s == 0 mod %s"%(self.polynomial(),
                                 mod_str)
        
    def _latex_(self):
        return latex(self.polynomial()) + " = " + latex(0) + \
            "\\text{ (mod }" + latex(self.modulus()) + "\\text{)}"

    def _cache_key(self):
        return 'CongruenceCondition', self.polynomial(), self.modulus()

class PowerCondition(PolynomialCondition):
    r"""
    Defines an expression that should be an n-th power.
    """

    def __init__(self, polynomial, least_exp=1):
        r"""
        Initializes a PowerCondition.

        INPUT:

        - ``polynomial`` -- A polynomial that
          should be an unknown power of some
          number.
        - ``least_exp`` -- A strictly positive integer
          (default: 1) that is the least power that
          this polynomial must be.
        """
        PolynomialCondition.__init__(self, polynomial)
        self._exp = least_exp

    def least_exponent(self):
        r"""
        OUTPUT:

        The smallest integer such that the polynomial in
        this condition is allowed to be that power of 
        some number.
        """
        return self._exp

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False,
                   verbose=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of those
        values for the variables such that the polynomial
        of this condition can be a power of some number.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.
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

        A pAdicTree object that contatins that part of the
        given pAdicTree which satisfies the polynomial of
        this condition being equal to some power, at least
        least_exponent, of some number.
        If complement is set to True will also give the
        complement of this tree in the given tree as a 
        second return value.
        """
        T1, T0 = PolynomialCondition.pAdic_tree(self, pAdic_tree=pAdic_tree,
                                                pAdics=pAdics, complement=True,
                                                verbose=verbose, precision=1,
                                                **kwds)
        Te = PolynomialCondition.pAdic_tree(self, pAdic_tree=pAdic_tree,
                                            pAdics=pAdics, complement=False,
                                            verbose=verbose, precision=self.least_exponent(),
                                            **kwds)
        T = T0.union(Te)
        if complement:
            return T, pAdic_tree.subtract(T)
        else:
            return T

    @cached_method
    def _x_str(self):
        i = 0
        x = 'x' + str(i)
        while x in self.variables():
            i += 1
            x = 'x' + str(i)
        return x

    @cached_method
    def _n_str(self):
        i = 0
        n = 'n' + str(i)
        while n in self.variables():
            i += 1
            n = 'n' + str(i)
        return n
        
    def _repr_(self):
        return "The condition that %s == %s^%s with %s >= %s"%(self.polynomial(),
                                                               self._x_str(),
                                                               self._n_str(),
                                                               self._n_str(),
                                                               self.least_exponent())

    def _repr_short(self):
        return "%s == %s^%s"%(self.polynomial(),
                              self._x_str(),
                              self._n_str())
        
    def _latex_(self):
        return latex(self.polynomial()) + " = " + \
            "x_{" + self._x_str()[1:] + "}" + \
            "^{n_{" + self._n_str()[1:] + "}}"

    def _cache_key(self):
        return 'PowerCondition', self.polynomial(), self.least_exponent()

class CoprimeCondition(Condition_base):
    r"""
    Defines the condition of variables being n-wise coprime.
    """

    def __init__(self, variables, n=2):
        r"""
        The constructor of a CoprimeCondition.

        INPUT:
        
        - ``variables`` -- A collection of variables on
          which this condition applies. This may be any
          form of a variable, but will be converted into
          strings. Multiple variables with the same name
          are therefore not very well supported and may
          cause unpredictable behavior.
        - ``n`` -- A non-negative integer (default: 2)
          indicating the size of subsets of the variables
          that should be coprime, e.g. n=2 means that
          the variables should be pairwise coprime and
          n=1 indicates all variables should be units.
        """
        Condition_base.__init__(self, variables)
        self._n = n

    def number_of_coprimes(self):
        r"""
        OUTPUT:

        An integer n such that this condition is that
        the variables should be n-wise coprime.
        """
        return self._n

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree such that
        all variables in this condition are n-wise coprime.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the
        given pAdicTree such that the variables in this
        condition are n-wise coprime.
        If complement was set to True will return a tuple
        with the afore mentioned as its first entry and the
        complement of that tree within the given pAdicTree
        as its second argument.
        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics,
                                   full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        T = pAdic_tree.root()
        tree_vars = pAdic_tree.variables()
        indices = tuple(tree_vars.index(var) for var in self.variables() if var in tree_vars)
        for node in T.children_at_level(1):
            if sum(c == 0 for i,c in enumerate(node.quotient_tuple()) if i in indices) >= self._n:
                node.remove()
        return pAdicTree(variables=pAdic_tree.variables(),
                         root=T)

    def _repr_(self):
        if self._n == 0:
            return "The condition that always holds"
        if self._n == 1:
            return "The condition that the variables %s are units."%(self.variables(),)
        if self._n == 2:
            return "The condition that the variables %s are pairwise coprime."%(self.variables(),)
        return "The condition that the variables %s are %s-wise coprime."%(self.variables(),
                                                                          self._n)

    def _repr_short(self):
        if self._n == 0:
            return "true"
        if self._n == 1:
            return "%s are units"%(self.variables(),)
        if self._n == 2:
            return "%s are pairwise coprime"%(self.variables(),)
        return "%s are %s-wise coprime"%(self.variables(),
                                          self._n)
        
    def _latex_(self):
        if self._n == 0:
            return "\\top"
        if self._n == 1:
            return latex(self.variables()) + "\\text{ are units}"
        if self._n == 2:
            return latex(self.variables()) + "\\text{ are pairwise coprime}"
        return latex(self.variables()) + "\\text{ are $" + str(self._n) + "$-wise coprime.}"

    def _cache_key(self):
        return 'CoprimeCondition', self.variables(), self._n

class NotCondition(Condition_base):
    r"""
    The condition that not another condition holds.
    """

    def __init__(self, other):
        r"""
        The constructor of a NotCondition.

        INPUT:

        - ``other`` -- The condition that this condition
          should be the negation of.
        """
        self._other = other
        Condition_base.__init__(self, other.variables())

    def negated_condition(self):
        r"""
        OUTPUT:

        The condition that this condition is the negation of.
        """
        return self._other

    def pAdic_tree(self, complement=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree that does not
        satisfy the condition stored in this condition.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the
        given pAdicTree which does not satisfy the condition
        defined by the condition stored in this condition.
        If complement was set to True will return a tuple with
        the afore mentioned as its first entry and the
        complement of that tree within the given pAdicTree as
        its second argument.
        """
        TY, TN = self._other.pAdic_tree(complement=True, **kwds)
        if complement:
            return TN, TY
        else:
            return TN

    def _repr_(self):
        s = self._other._repr_()
        s = s.replace(' are ', ' are not ')
        s = s.replace(' is ', ' is not ')
        s = s.replace(' not not ', ' ')
        s = s.replace(' always ', '<tmp>')
        s = s.replace(' never ', ' always ')
        s = s.replace('<tmp>', ' never ')
        s = s.replace('==', '<tmp>')
        s = s.replace('~=', '==')
        s = s.replace('<tmp>', '~=')
        s = s.replace(' and ', '<tmp>')
        s = s.replace(' or ', ' and ')
        s = s.replace('<tmp>', ' or ')
        return s

    def _repr_short(self):
        s = self._other._repr_short()
        s = s.replace(' are ', ' are not ')
        s = s.replace(' is ', ' is not ')
        s = s.replace(' not not ', ' ')
        s = s.replace(' always ', '<tmp>')
        s = s.replace(' never ', ' always ')
        s = s.replace('<tmp>', ' never ')
        s = s.replace('true', '<tmp>')
        s = s.replace('false', 'true')
        s = s.replace('<tmp>', 'false')
        s = s.replace('==', '<tmp>')
        s = s.replace('~=', '==')
        s = s.replace('<tmp>', '~=')
        s = s.replace(' and ', '<tmp>')
        s = s.replace(' or ', ' and ')
        s = s.replace('<tmp>', ' or ')
        return s
        
    def _latex_(self):
        s = self._other._latex_()
        s = s.replace(' are ', ' are not ')
        s = s.replace(' is ', ' is not ')
        s = s.replace(' not not ', ' ')
        s = s.replace('=', '<tmp>')
        s = s.replace('\\neq', '=')
        s = s.replace('<tmp>', '\\neq')
        s = s.replace('\\wedge', '<tmp>')
        s = s.replace('\\vee', '\\wedge')
        s = s.replace('<tmp>', '\\vee')
        s = s.replace('\\top', '<tmp>')
        s = s.replace('\\bot', '\\top')
        s = s.replace('<tmp>', '\\bot')
        return s

    def _cache_key(self):
        return 'NotCondition', self._other

class AndCondition(Condition_base):
    r"""
    The condition that both conditions hold.
    """

    def __init__(self, left, right):
        r"""
        The constructor of an AndCondition.

        INPUT:

        - ``left`` -- The first condition that should
          hold for this condition to hold.
        - ``right`` - The second condition that should
          hold for this condition to hold.
        """
        self._left = left
        self._right = right
        variables = list(left.variables())
        for var in right.variables():
            if var not in variables:
                variables.append(var)
        Condition_base.__init__(self, variables)

    def other_condition(self):
        r"""
        OUTPUT:

        The conditions that should hold for this condition
        to hold.
        """
        return self._left, self._right

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree that satisfies
        both conditions stored in this condition.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the
        given pAdicTree which satisfies both conditions
        defined in this condition.
        If complement was set to True will return a tuple with
        the afore mentioned as its first entry and the
        complement of that tree within the given pAdicTree as
        its second argument.
        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics,
                                   full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        T1 = self._left.pAdic_tree(pAdic_tree=pAdic_tree,
                                   pAdics=pAdics, complement=False,
                                   **kwds)
        T2 = self._right.pAdic_tree(pAdic_tree=pAdic_tree,
                                    pAdics=pAdics, complement=False,
                                    **kwds)
        T = T1.intersection(T2)
        if complement:
            return T, pAdic_tree.difference(T)
        else:
            return T

    def _repr_(self):
        right_str = self._right._repr_()
        right_str = right_str[0].lower() + right_str[1:]
        left_str = self._left._repr_()
        if left_str.endswith('.'):
            left_str = left_str[:-1]
        return left_str + " and " + right_str

    def _repr_short(self):
        return self._left._repr_short() + " and " + self._right._repr_short()
        
    def _latex_(self):
        return self._left._latex_() + " \\wedge " + self._right._latex_()

    def _cache_key(self):
        return 'AndCondition', self._left, self._right

class OrCondition(Condition_base):
    r"""
    The condition that either one of two conditions holds.
    """

    def __init__(self, left, right):
        r"""
        The constructor of an OrCondition.

        INPUT:

        - ``left`` -- The first condition that could
          hold for this condition to hold.
        - ``right`` - The second condition that could
          hold for this condition to hold.
        """
        self._left = left
        self._right = right
        variables = list(left.variables())
        for var in right.variables():
            if var not in variables:
                variables.append(var)
        Condition_base.__init__(self, variables)

    def other_condition(self):
        r"""
        OUTPUT:

        The conditions that could hold for this condition
        to hold.
        """
        return self._left, self._right

    def pAdic_tree(self, complement=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.
        
        Given a p-adic tree, returns the subtree of which
        each value satisfy one of the two conditions in this
        condition.
        
        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.

        OUTPUT:

        A pAdicTree object that contains that part of the
        given pAdicTree of which the value satisfy at least
        one of the two conditions defined in this condition.
        If complement was set to True will return a tuple with
        the afore mentioned as its first entry and the
        complement of that tree within the given pAdicTree as
        its second argument.
        """
        if pAdic_tree is None:
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics,
                                   full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        T1 = self._left.pAdic_tree(pAdic_tree=pAdic_tree,
                                   pAdics=pAdics, complement=False,
                                   **kwds)
        T2 = self._right.pAdic_tree(pAdic_tree=pAdic_tree,
                                   pAdics=pAdics, complement=False,
                                   **kwds)
        T = T1.union(T2)
        if complement:
            return T, pAdic_tree.difference(T)
        else:
            return T

    def _repr_(self):
        right_str = self._right._repr_()
        right_str = right_str[0].lower() + right_str[1:]
        left_str = self._left._repr_()
        if left_str.endswith('.'):
            left_str = left_str[:-1]
        return left_str + " or " + right_str

    def _repr_short(self):
        return self._left._repr_short() + " or " + self._right._repr_short()
        
    def _latex_(self):
        return self._left._latex_() + " \\vee " + self._right._latex_()

    def _cache_key(self):
        return 'OrCondition', self._left, self._right

class TreeCondition(Condition_base):
    r"""
    A condition that the value should be part of some pAdicTree.
    """

    def __init__(self, pAdic_tree):
        r"""
        Te constructor of a TreeCondition.

        INPUT:

        - ``pAdic_tree`` -- A pAdicTree that contains the
          variables and values they should attain for this
          condition to hold.
        """
        self._T = pAdic_tree
        Condition_base.__init__(self, pAdic_tree.variables())

    def pAdic_tree(self, pAdic_tree=None, pAdics=None, complement=False, **kwds):
        r"""
        Gives this condition as a pAdicTree.

        INPUT:
        
        - ``pAdic_tree`` -- A pAdicTree object (default:None)
          on which this condition should be applied. If set
          to None will be initiated as the full tree with
          the given pAdics.
        - ``pAdics`` -- A pAdicBase object (default: None)
          determining the pAdics that should be used. If set
          to None will use the pAdics of the given pAdicTree
          instead. If that is also set to None, will use the
          pAdics of the tree stored in this Condition instead.
        - ``complement`` -- A boolean (default: False)
          determining whether the complement of the result
          should be returned.

        OUTPUT:

        If the given pAdicTree has no common pAdics with the
        pAdicTree stored in this Condition will return the given
        pAdicTree. Otherwise will return a pAdicTree with the
        same variables as the given pAdicTree and values such
        that for each value there is a value in the pAdicTree
        stored in this Condition with the same value for each
        (common) variable.

        If complement was set to True will return a tuple with
        the afore mentioned as its first entry and the
        complement of that tree within the given pAdicTree as
        its second argument.
        """
        if pAdic_tree is None:
            if pAdics is None:
                if complement:
                    return self._T, self._T.complement()
                else:
                    return self._T
            pAdic_tree = pAdicTree(variables=self.variables(),
                                   pAdics=pAdics,
                                   full=True)
        if pAdics is None:
            pAdics = pAdic_tree.pAdics()
        if pAdics == self._T.pAdics():
            result = self._T.intersection(pAdic_tree).change_variables_to(pAdic_tree.variables())
            if complement:
                return result, pAdic_tree.difference(result)
            else:
                return result
        elif complement:
            return pAdic_tree, pAdic_tree
        else:
            return pAdic_tree

    def _repr_len(self, max_item=50, max_char=1000):
        r"""
        Gives a string representation of this Condition of at
        most a given length

        INPUT:

        - ``max_item`` -- A non-negative integer (default: 10)
          indicating the maximal number of items to be included
          in this representation.
        - ``max_char`` -- A non-negative integer (default: 200)
          giving the maximal number of character to be used in
          the string representation of this object.
        """
        l = len(self.variables())
        if l == 0:
            return "true"
        values, modulus = self._T.give_as_congruence_condition()
        if hasattr(modulus, 'is_principal') and modulus.is_principal():
            modulus = modulus.gens_reduced()[0]
        if len(values) <= max_item:
            result = str(self.variables() if l > 1 else self.variables()[0]) + \
                     " == "
            for i, value in enumerate(values):
                if len(result) > max_char:
                    break
                if i > 0:
                    result += ", "
                result += str(value if l > 1 else value[0])
            result += " mod " + \
                  (modulus._repr_short() if hasattr(modulus, '_repr_short') else str(modulus))
            if len(result) <= max_char:
                return result
        return str(self.variables() if l > 1 else self.variables()[0]) + \
               " is 1 of " + \
               str(len(values)) + \
               " possibilities mod " + \
               (modulus._repr_short() if hasattr(modulus, '_repr_short') else str(modulus))
        
    def _repr_(self):
        result = "The condition that " + self._repr_len()
        result.replace("true", "always holds")
        result.replace("mod", "modulo")
        return result
               
    def _repr_short(self):
        return self._repr_len(max_item=10, max_char=40)

    def _latex_(self):
        l = len(self.variables())
        if l == 0:
            return "\\top"
        values, modulus = self._T.give_as_congruence_condition
        if hasattr(modulus, 'is_principal') and modulus.is_principal():
            modulus = modulus.gens_reduced()[0]
        result = latex(self.variables() if l > 1 else self.variables()[0]) + \
                 " = "
        for i, value in enumerate(values):
            if i > 0:
                result += ", "
            result += latex(value if l > 1 else value[0])
        result += " \\text{ (mod }" + \
                  latex(modulus) + \
                  "\\text{)}"
        return result

    def _cache_key(self):
        return 'TreeCondition', self._T

class ConditionalValue(SageObject):
    r"""
    Some value that depends on some condition.
    """

    def __init__(self, val_con):
        r"""
        Initializes this object.

        INPUT:

        - ``val_con`` - A list of tuples, where each
          tuple consists of a value and a condition
          on when this value is attained, in that order.
          A value can be any object, whilst a condition
          must extend Condition_base. The different
          conditions do not have to include all
          possibilities, nor do they have to be exclusive
          of one another, but not adhering to this will
          reflect in the resulting object.
        """
        self._vals = tuple(vc[0] for vc in val_con)
        self._con = tuple(vc[1] for vc in val_con)
        for c in self._con:
            if not isinstance(c, Condition_base):
                raise ValueError("%s is not a condition"%(c,))

    def no_value(self):
        r"""
        Tells if this conditional value does not attain any value.

        OUTPUT:

        True - If this condition does not contain any possible value.
        False - If this condition contains at least one value.
        """
        return len(self._vals) == 0

    def _repr_lines(self):
        if self.no_value():
            return []
        result = [str(val) for val in self._vals]
        l = max(len(r) for r in result) + 1
        result = [r +
                  (' '*(l - len(r))) +
                  "if " +
                  self._con[i]._repr_short()
                  for i, r in enumerate(result)]
        return result
    
    def _repr_(self):
        if self.no_value():
            return "Conditional value with no possible value."
        lines = self._repr_lines()
        result = ""
        for i, line in enumerate(lines):
            if i > 0:
                result += "\n"
            result += line
        return result

    def _latex_lines(self):
        if self.no_value():
            return []
        return [latex(val) +
                "& \\text{ if }" +
                latex(self._con[i])
                for i, val in enumerate(self._vals)]
    
    def _latex_(self):
        if self.no_value():
            return "\\text{Conditional value with no possible value.}"
        result = "\\left\\{ \\begin{array}{lr}\n"
        for i, line in enumerate(self._latex_lines()):
            if i > 0:
                result += " \\\\\n"
            result += line
        result += "\n \\end{array} \\right."
        return result

    def __add__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR, self, other)

    def __radd__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR, other, self)

    def __sub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, self, other)

    def __rsub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, other, self)

    def __mul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, self, other)

    def __rmul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, other, self)

    def __div__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, self, other)

    def __rdiv__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, other, self)

    def __pow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, self, other)

    def __rpow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, other, self)

    def __iter__(self):
        return iter(zip(self._vals, self._con))

    def __len__(self):
        return len(self._vals)

    def __getitem__(self, index):
        return (self._vals[index], self._con[index])

    def _cache_key(self):
        return tuple((val, 'if', con) for val, con in self)

class ConditionalExpression(SageObject):
    SUM_OPERATOR = ('+', '+', 0)
    MINUS_OPERATOR = ('-', '-', 0.5)
    PRODUCT_OPERATOR = ('*', '\cdot', 2)
    DIVISION_OPERATOR = ('/','/', 2.5)
    EXPONENT_OPERATOR = ('^', '^', 4.5)
    r"""
    An expression containing conditional values.
    """
    
    def __init__(self, operator, left, right):
        r"""
        Initializes a conditional expression.

        INPUT:

        - ``operator`` -- A tuple containing in this order
          a string representing the operator, a string
          that will produce the operator in latex and a
          non-negative integer indicating the power of the
          operator.
        """
        self._op = operator
        self._left = left
        self._right = right

    def left(self):
        r"""
        Gives the left side of this expression.

        A ConditionalExpression consists of a left and
        a right side seperated by some operator.

        OUTPUT:

        The left side of this expression.
        """
        return self._left

    def operator(self):
        r"""
        Gives the operator of this expression.

        A ConditionalExpression consists of a left and
        a right side seperated by some operator.

        OUTPUT:

        The left side of this expression.
        """
        return self._op

    def right(self):
        r"""
        Gives the right side of this expression.

        A ConditionalExpression consists of a left and
        a right side seperated by some operator.

        OUTPUT:

        The right side of this expression.
        """
        return self._right

    def _factor_side(self, side):
        r"""
        Factors a side. See meth::factor for details.
        """
        if isinstance(side, ConditionalExpression):
            return side.factors()
        if hasattr(side, 'factor'):
            F = side.factor()
            result = {f: e for f,e in F}
            if hasattr(F, 'unit') and F.unit() != 1:
                result[F.unit()] = 1
            return result
        return {side: 1}
    
    def factors(self):
        r"""
        Gives the factors in this expression and their exponents.

        The factoring will work as follows:
        - If the operator is a PRODUCT_OPERATOR will factor left
          and right and combine both by taking the union of the
          factors of left and right and taking as the exponent
          the sum of the exponents if the factor appeared at both
          sides or the exponent of the side it appeared on otherwise.
        - If the operator is a DIVISION_OPERATOR will factor left
          and right and combine both by taking the union of the
          factors of left and right and taking as exponents the
          exponent of that factor in left minus the exponent of
          that factor in right. If the factor appeared only on
          the left will take the corresponding exponent instead
          and if the factor only appeared on the right will take
          0 minus the corresponding exponent instead.
        - If the operator is a EXPONENT_OPERATOR will factor left
          and take all those factors. The exponent of each factor
          will be the corresponding exponent in left times the
          right side.
        - If the operator is anything else, will return this as
          the only factor and 1 as the exponent.
        - When factoring a side, will use the method factors if
          the side is a ConditionalExponent, will attempt to
          use the method factor if it is any other object with
          such a method and in any other case will assume that
          side to be the factor with exponent 1. If the method
          factor is used and the resulting object has a method
          unit, the result of this method will be added as a
          factor with respective exponent 1 unless it is 1
          itself.

        OUTPUT:

        A dictionary containing as keys the factors of this
        expression and as values their corresponding exponents.
        Note that these factors and exponents could be any
        type of expression that could normally appear in a
        conditional expression.
        """
        if self._op == ConditionalExpression.PRODUCT_OPERATOR:
            result = self._factor_side(self._left)
            extra = self._factor_side(self._right)
            for f in extra:
                if f in result and result[f] != 0:
                    result[f] = result[f] + extra[f]
                elif extra[f] != 0:
                    result[f] = extra[f]
            return result
        if self._op == ConditionalExpression.DIVISION_OPERATOR:
            result = self._factor_side(self._left)
            extra = self._factor_side(self._right)
            for f in extra:
                if f in result and result[f] != 0:
                    result[f] = result[f] - extra[f]
                elif extra[f] != 0:
                    result[f] = 0 - extra[f]
            return result
        if self._op == ConditionalExpression.EXPONENT_OPERATOR:
            result = self._factor_side(self._left)
            for f in result:
                if result[f] == 1:
                    result[f] = self._right
                elif self._right != 1:
                    result[f] = result[f] * self._right
            return result
        return {self: 1}

    def _repr_side(self, side, vals, bracket_level):
        if isinstance(side, ConditionalExpression):
            return side._repr_info(vals, bracket_level)
        if isinstance(side, ConditionalValue):
            vals.append(side)
            return "n" + str(len(vals) - 1)
        if hasattr(side, "_repr_short"):
            return side._repr_short()
        return str(side)

    def _operation_on(self, left, right):
        r"""
        Performs the operation of this expression on two values.

        Depending on the operation will do the following
        - If the operation is SUM_OPERATION will return the sum
          of left plus right
        - If the operation is MINUS_OPERATION will return left
          minus right.
        - If the operation is PRODUCT_OPERATION will return the
          product of left with right.
        - If the operation is DIVISION_OPERATION will return
          left divided by right.
        - If the operation is EXPONENT_OPERATION will return
          left to the power right
        - If the operation is anything else will return a
          ValueError.
        """
        if self._op == ConditionalExpression.SUM_OPERATOR:
            return left + right
        if self._op == ConditionalExpression.MINUS_OPERATOR:
            return left - right
        if self._op == ConditionalExpression.PRODUCT_OPERATOR:
            return left * right
        if self._op == ConditionalExpression.DIVISION_OPERATOR:
            return left / right
        if self._op == ConditionalExpression.EXPONENT_OPERATOR:
            return left ^ right
        raise ValueError('Can not evaluate operation %s'%(self._op,))
    
    def value(self):
        r"""
        Turns this conditional expression into a conditional value.

        Attempts to turn this conditional expression into a conditional
        value by taking for each conditional value in the expression
        a list of possible combinations of their conditions as the
        conditions for the resulting conditional value. For each of
        these conditions will substitute the appropiate value
        for each conditional value in this expression and evaluate
        the expression assigning the result as the corresponding
        value of the resulting conditional value.

        OUTPUT:

        A ConditionalValue or some element that has exactly the same
        values as this expression on every choice of parameters.
        """
        left = self._left
        if isinstance(left, ConditionalExpression):
            left = left.value()
        right = self._right
        if isinstance(right, ConditionalExpression):
            right = right.value()
        if isinstance(left, ConditionalValue):
            result = []
            for val_l, con_l in left:
                if isinstance(right, ConditionalValue):
                    for val_r, con_r in right:
                        result.append((self._operation_on(val_l, val_r),
                                       con_l & con_r))
                else:
                    result.append((self._operation_on(val_l, right),
                                   con_l))
            return ConditionalValue(result)
        if isinstance(right, ConditionalValue):
            result = []
            for val, con in right:
                result.append((self._operation_on(left, val),
                               con))
            return ConditionalValue(result)
        return self._operation_on(left, right)
        
    def _repr_info(self, vals, bracket_level):
        result = ""
        if bracket_level > self._op[2]:
            result += "("
        result += self._repr_side(self._left, vals, floor(self._op[2]))
        result += self._op[0]
        result += self._repr_side(self._right, vals, ceil(self._op[2]))
        if bracket_level > self._op[2]:
            result += ")"
        return result

    def _repr_(self):
        vals = []
        result = self._repr_info(vals, 0)
        if len(vals) > 0:
            result += "\n where \n"
            front_space = ceil(ZZ(len(vals)).log(10)) + 5
            for i, val in enumerate(vals):
                lines = val._repr_lines()
                if len(lines) == 0:
                    if i > 0:
                        result += "\n"
                    result += "n" + str(i) + " ="
                    result += " no possible value"
                for j, line in enumerate(lines):
                    if i + j > 0:
                        result += "\n"
                    r = ""
                    if j == 0:
                        r += "n" + str(i) + " ="
                    result += r + " "*(front_space-len(r)) + \
                              line
        return result

    def _latex_side(self, side, vals, bracket_level):
        if isinstance(side, ConditionalExpression):
            return side._latex_info(vals, bracket_level)
        if isinstance(side, ConditionalValue):
            vals.append(side)
            return "n_{" + str(len(vals) - 1) + "}"
        return str(side)
        
    def _latex_info(self, vals, bracket_level):
        result = ""
        if bracket_level > self._op[2]:
            result += "\\left("
        result += self._latex_side(self._left, vals, floor(self._op[2]))
        result += self._op[1]
        result += self._latex_side(self._right, vals, ceil(self._op[2]))
        if bracket_level > self._op[2]:
            result += "\\right)"
        return result

    def _latex_(self):
        vals = []
        result = self._latex_info(vals,0)
        result += "\n\\\\ \text{ where } \\\\\n"
        for i, val in enumerate(vals):
            if i > 0:
                result += " \\\n"
            result += "n_{" + str(i+1) + "} = "
            result += " & " + latex(val)
        return result

    def __add__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR, self, other)

    def __radd__(self, other):
        return ConditionalExpression(ConditionalExpression.SUM_OPERATOR, other, self)

    def __sub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, self, other)

    def __rsub__(self, other):
        return ConditionalExpression(ConditionalExpression.MINUS_OPERATOR, other, self)

    def __mul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, self, other)

    def __rmul__(self, other):
        return ConditionalExpression(ConditionalExpression.PRODUCT_OPERATOR, other, self)

    def __div__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, self, other)

    def __rdiv__(self, other):
        return ConditionalExpression(ConditionalExpression.DIVISION_OPERATOR, other, self)

    def __pow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, self, other)

    def __rpow__(self, other):
        return ConditionalExpression(ConditionalExpression.EXPONENT_OPERATOR, other, self)

def apply_to_conditional_value(function, value, singleton=False,
                               use_condition=False, default_condition=None):
    r"""
    Applies a function to a conditional value.

    INPUT:
    
    - ``function`` -- Any function
    - ``value`` -- Any value that is accepted by the
      given function as an input or a ConditionalValue
      that contains such values.
    - ``singleton`` -- A boolean value (default: False)
      indicating whether a ConditionalValue with only
      one possibility should be returned.
    - ``use_condition`` -- A boolean value (default:
      False). If set to true, will pass a value
      condition pair to the function for each possible
      value instead of only the value. This is useful
      for functions which also require the condition
      on which this value depends.
    - ``default_condition`` -- If the given value was
      not ConditionalValue, will use this to determine
      the condition corresponding to the value if it
      is required to pass to be passed to the funcion
      or to construct a ConditionalValue at the end.

    OUTPUT:

    The function evaluated on the given value. If the
    given value was a ConditionalValue this means the function
    is evaluated on every value in the ConditionalValue
    producing a new ConditionalValue of possible outcomes.
    Conditions which produce the same outcome will be combined
    using the OrCondition. If the resulting ConditionalValue
    would have only one possibility and singleton is set to
    False, will return the value of that single possibility
    instead of the whole ConditionalValue.
    """
    if isinstance(value, ConditionalValue):
        result = {}
        for val, con in value:
            if use_condition:
                f_val = function(val, con)
            else:
                f_val = function(val)
            if f_val in result:
                result[f_val] = (result[f_val] | con)
            else:
                result[f_val] = con
        if not singleton and len(result) == 1:
            return result.keys()[0]
        else:
            return ConditionalValue(list(result.iteritems()))
    else:
        if use_condition:
            result = function(value, default_condition)
        else:
            result = function(value)
        if singleton and default_conditon != None:
            return ConditionalValue([(result, default_condition)])
        else:
            return result

def conditional_product(*args):
    r"""
    Given multiple ConditionalValues creates a
    ConditionalValue of lists of possible values.

    INPUT:
    
    Any amount of arguments, all of which can be
    an instance of ConditionalValue or any other
    value.

    OUTPUT:

    A ConditionalValue of which the values are all
    possible lists with as i-th entry a possible
    value of the i-th given ConditionalValue. For
    each value the corresponding condition is the
    condition that all the corresponding conditions
    for each entry in the respective ConditionalValue
    hold.

    Note that if each entry given is not a
    ConditionalValue, the result will also not be
    a ConditionalValue.
    """
    if len(args) == 0:
        raise ValueError("conditional_zip requires at least one argument.")
    for i in range(len(args)):
        if not isinstance(args[i], ConditionalValue):
            args[i] = [(args[i], None)]
    # Helper function
    def combine_conditions(C1, C2):
        if C1 is None:
            return C2
        if C2 is None:
            return C1
        return C1 & C2
    # Make a list of lists of value condition pairs
    result = itertools.product(*args)
    # Turn into a list of pairs of a list of values
    # and the corresponding list of conditions.
    result = [zip(*val_con) for val_con in result]
    # Turn into a list of pairs of a list of values
    # and the corresponding condition, combined with and.
    result = [(val, reduce(combine_conditions, con, None)) for val, con in result]
    if len(result) == 1:
        return result[0]
    else:
        return ConditionalValue(result)
