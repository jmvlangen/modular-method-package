r"""
Basic methods for dealing with p-adic completions of number fields

This file contains the class pAdicBase which is a wrapper for easy
access to functionality needed by more advanced p-adic computational
code.

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Joey van Langen (2018-07-13): initial version

"""

# ****************************************************************************
#       Copyright (C) 2018 Joey van Langen <j.m.van.langen@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.number_field.order import Order
from sage.rings.number_field.number_field import NumberField_generic
from sage.structure.sage_object import SageObject

class pAdicBase(SageObject):
    r"""
    The class ``pAdicBase`` stores basic information that forms
    the basis of p-adic arithmatic. It contains methods to retrieve
    the following data:
        - a number field over which we are working p-adicly
        - the ring of integers of this number field
        - a prime ideal of this ring for which we are working p-adicly
        - a uniformizer for the prime ideal
        - the residue field of this prime ideal
        - the characteristic of the residue field
        - representatives for the residue field
        - the quotient ring of the ring of integers modulo some power of the
          prime ideal
        - the valuation provided by the prime ideal
        - the zero element in this context.
        
    EXAMPLES::
        
        sage: pAdics = pAdicBase(ZZ, 5)
        sage: pAdics
        p-adic information of Rational Field with respect to (5).
        sage: pAdics.valuation(100)
        2
        
    Example using a number field::
    
        sage: K = QuadraticField(-7)
        sage: pAdics = pAdicBase(K, K.prime_above(2))
        sage: pAdics
        p-adic information of Number Field in a with defining polynomial x^2 + 7 with respect to (-1/2*a + 1/2).
        sage: pAdics.characteristic()
        2
    """
    
    def __init__(self, ring, prime):
        r"""
        Initialization of the ``pAdicBase`` class.
        
        INPUT:
        
        - ``ring`` -- A number field or order therein.
        - ``prime`` -- A prime ideal of the maximal order of the
          number field given by the first argument or a generator
          thereof.
          
        EXAMPLES:
            
            sage: pAdicBase(QQ,3)
            p-adic information of Rational Field with respect to (3).
            
        Example using a number field::
        
            sage: K = CyclotomicField(5)
            sage: pAdicBase(K, K.prime_above(5))
            p-adic information of Cyclotomic Field of order 5 and degree 4 with respect to (zeta5 - 1).
            
        """
        if isinstance(ring, Order):
            ring = ring.number_field()
        if ring == ZZ:
            ring = QQ
        if isinstance(ring, NumberField_generic):
            self._R = ring.maximal_order()
        elif ring == QQ:
            self._R = ZZ
        else:
            raise ValueError("%s is not a number field."%ring)
        if prime in self._R:
            prime = self._R.ideal(prime)
        if prime not in self._R.ideal_monoid() or not prime.is_prime():
            raise ValueError("%s is not a prime ideal of %s."%(prime, self._R))
        self._P = prime
        self._ext = dict() #Residue field as vector space, see extension_vector_space
        
    def number_field(self):
        r"""
        Returns the number field associated to these p-adics.
        
        OUTPUT:
        The number field associated to these p-adics.
        
        EXAMPLE::
        
            sage: K = QuadraticField(7)
            sage: pAdics = pAdicBase(K, K.prime_above(2))
            sage: K == pAdics.number_field()
            True
            
        """
        if self._R == ZZ:
            return QQ
        else:
            return self._R.number_field()
    
    def order(self):
        r"""
        Returns the maximal order associated to these p-adics.
        
        OUTPUT:
        The maximal order in the number field returned by
        :func:`number_field`.
        
        EXAMPLES::
        
            sage: K = QuadraticField(-5)
            sage: pAdicBase(K, K.prime_above(17)).order()
            Maximal Order in Number Field in a with defining polynomial x^2 + 5
            
        ::
        
            sage: K = CyclotomicField(7)
            sage: pAdics = pAdicBase(K, K.prime_above(3))
            sage: pAdics.order() == K.ring_of_integers()
            True
            
        """
        return self._R
        
    def prime_ideal(self):
        r"""
        Returns the prime ideal associated to these p-adics
        
        OUTPUT:
        The prime ideal associated to these p-adics.
        
        EXAMPLES::
        
            sage: pAdics = pAdicBase(ZZ,13)
            sage: pAdics.prime_ideal()
            Principal ideal (13) of Integer Ring
            
        ::
        
            sage: K = QuadraticField(-1)
            sage: pAdics = pAdicBase(K, K.prime_above(7))
            sage: pAdics.prime_ideal() == K.prime_above(7)
            True
            
        """
        return self._P

    def prime_below(self, R):
        r"""
        Gives the prime that lies below the prime of this
        pAdicBase.

        INPUT:
        
        - ``R`` -- A subring of the number field of this
          pAdicBase.

        OUTPUT:

        A prime of the ring of fractions of R that lies
        below the prime stored in this pAdicBase.
        """
        K = R.fraction_field()
        L = self.number_field()
        Q = self.prime_ideal()
        if L is QQ:
            p = Q.gens()[0]
        else:
            p = Q.smallest_integer()
        if K == QQ:
            return p
        ls = K.gen().minpoly().change_ring(L).roots()
        if len(ls) < 1:
            raise ValueError("%s is not a subring of %s"%(R,
                                                          L))
        iota = K.hom([ls[0][0]], L)
        for P in K.primes_above(p):
            PL = L.ideal([iota(g) for g in P.gens()])
            if Q in L.primes_above(PL):
                return P
        raise ValueError("No prime in %s lies below %s"%(K, Q))
        
    def uniformizer(self):
        r"""
        Returns a uniformizer associated to these p-adics.
        
        OUTPUT:
        An element of the prime ideal returned by
        :func:`prime_ideal` that has valuation 1 with
        respect to this prime ideal.
        
        EXAMPLE:
        
            sage: K.<a> = QuadraticField(10)
            sage: P = K.ideal(3, a + 1)
            sage: pAdics = pAdicBase(K, P)
            sage: pi = pAdics.uniformizer(); pi
            3
            sage: pAdics.valuation(pi)
            1
        """
        if not hasattr(self," _pi"):
            if self._R == ZZ:
                self._pi = self._P.gens()[0]
            else:
                for g in self._P.gens():
                    if g.ord(self._P) == 1:
                        self._pi = g
        return self._pi
            
    def valuation(self, x):
        r"""
        Calculates the valuation of some element.
        
        INPUT:
        
        - ``x`` -- An element of the number field returned by
          :func:`number_field`.
        
        OUTPUT:
        The valuation of ``x`` with respect tot the valuation
        associated to these p-adics, i.e. the order of the
        prime ideal returned by :func:`prime_ideal` in the
        prime factorization of ``x``.
        
        EXAMPLES::
        
            sage: pAdics = pAdicBase(QQ, 2)
            sage: pAdics.valuation(24)
            3
            sage: pAdics.valuation(3/2)
            -1
            
        ::
            
            sage: K.<a> = QuadraticField(-3)
            sage: pAdics = pAdicBase(K, K.prime_above(7))
            sage: pAdics.valuation(18*a + 20)
            3
            sage: pAdics.valuation((9*a - 3)/14)
            -1
            
        """
        x = self.number_field()(x)
        if self._R == ZZ:
            return x.ord(self.uniformizer())
        else:
            return x.ord(self.prime_ideal())
            
    def _power_series(self, x, prec, F, pi):
        if prec > 0:
            c0 = F.lift(F(x))
            return [c0] + self._power_series((x-c0)/pi, prec-1, F, pi)
        else:
            return []
            
    def power_series(self, x, precision=10):
        """r
        Writes an element x as a power series in the uniformizer.
        
        INPUT:
        
        - ``x`` -- An element of the number field returned
          by :func:`number_field`, which has non-negative
          valuation with respect to the prime ideal returned
          by :func:`prime_ideal`.
        - ``precision`` -- A strictly positive integer
          (default: 10) which indicates to what order the
          power series should extend. Note that precision
          10 will result in a power series with highest
          power 
          
        OUTPUT:
        A list of elements `c_i`, such that
        - Each `c_i` is a representative as returned by
          the function :func:`representatives`.
        - We have that
        .. MATH::
        
            \sum_{i=0}^{k-1} c_i * \pi^i \equiv x \text{(mod} P^k),
            
        where `P` is the prime ideal returned by :func:`prime_ideal`
        and `k` is the given argument precision.
        """
        if self.valuation(x) < 0:
            raise ValueError("%s is not integral with respect to %s."%(x,
                                                                       self.prime_ideal()))
        if precision not in ZZ or precision <= 0:
            raise ValueError("%s is not a strictly positive integer."%(precision))
        return self._power_series(x, precision, self.residue_field(),
                                  self.uniformizer())
        
    def residue_field(self):
        r"""
        Returns the residue field associated to these p-adics.
        
        OUTPUT:
        The quotient of the order returned by :func:`order`
        modulo the prime returned by :func:`prime_ideal`,
        as a residue field object.
        
        EXAMPLES::
        
            sage: pAdicBase(ZZ, 5).residue_field()
            Residue field of Integers modulo 5
            
        ::
        
            sage: K = QuadraticField(-1)
            sage: pAdicBase(K, K.prime_above(5)).residue_field()
            Residue field of Fractional ideal (-a - 2)
            
        """
        if not hasattr(self, "_F"):
            self._F = self.number_field().residue_field(self._P)
        return self._F
    
    def representatives(self, width=1):
        r"""
        Returns a generator of representatives for a module of
        the residue field associated to these p-adics.
        
        INPUT:
        
        - ``width`` -- A strictly positive integer (default: `1`)
          that determines the rank of the module.
          
        OUTPUT:
        A generator of tuples of length ``width`` containing elements
        of the order returned by :func:`order`, which form
        representatives of F^width, where F is the residue field
        returned by :func:`residue_field`
        
        EXAMPLES:
        
        Determining the squares modulo 11. Note that the function
        returns tuples, hence we have to take the first argument::
        
            sage: pAdics = pAdicBase(ZZ, 11)
            sage: F = pAdics.residue_field()
            sage: squares = []
            sage: for x in pAdics.representatives():
            ....:     x2 = F(x[0]^2)
            ....:     if x2 not in squares:
            ....:         squares.append(x2)
            ....: squares.sort();
            ....: print squares
            ....:
            [1, 3, 4, 5, 9]
            
        Calculating representatives of possible roots modulo a
        prime ideal for a polynomial in two variables::
        
            sage: K = QuadraticField(5)
            sage: pAdics = pAdicBase(K, K.prime_above(11))
            sage: R = pAdics.order()
            sage: S.<x,y> = R[]
            sage: f = x^3 - y^2
            sage: roots = []
            sage: for z in pAdics.representatives(width=2):
            ....:     if pAdics.valuation(f(z)) > 0:
            ....:         roots.append(z)
            ....: print roots
            ....:
            [(0, 0), (7/2*a + 7/2, 7/2*a + 7/2), (a + 1, 3/2*a + 3/2), (3*a + 3, 5*a + 5), (5*a + 5, 3*a + 3), (4*a + 4, a + 1), (4*a + 4, 9/2*a + 9/2), (5*a + 5, 5/2*a + 5/2), (3*a + 3, 1/2*a + 1/2), (a + 1, 4*a + 4), (7/2*a + 7/2, 2*a + 2)]
            
        """
        if width not in ZZ or width <= 0:
            raise ValueError("width should be a strictly positive integer, not %s."%(width,))
        for cfs in self.residue_field()**width:
            yield tuple([self.residue_field().lift(c) for c in cfs])
            
    def size_residue_field(self):
        r"""
        Returns the size of the residue field associated to
        these p-adics.
        
        OUTPUT:
        The cardinality of the residue field returned by
        :func:`residue_field`.
        
        EXAMPLE::
        
            sage: K = QuadraticField(3)
            sage: pAdics = pAdicBase(K, K.prime_above(5))
            sage: pAdics.size_residue_field()
            25
            
        """
        return self.residue_field().cardinality()
        
    def characteristic(self):
        r"""
        Returns the characteristic of the residue field
        associated to these p-adics.
        
        OUTPUT:
        The characteristic of the residue field returned by the
        function :func:`residue_field`, i.e. the prime number in
        ZZ above the prime returned by :func:`prime_ideal`
        
        EXAMPLE::
            
            sage: K = CyclotomicField(11)
            sage: pAdics = pAdicBase(K, K.prime_above(7))
            sage: pAdics.characteristic()
            7
        """
        return self.residue_field().characteristic()
        
    def quotient_ring(self, n, names='param'):
        r"""
        Returns the quotient ring by a power of the prime ideal
        associated to these p-adics.
        
        INPUT:
        
        - ``n`` -- A non-negative integer determining which power of
          the prime ideal should be quotiented out.
        - ``names`` -- An identifier or list of identifiers 
          (default: 'param') that can be used as the names of quotient
          variables.
          
        OUTPUT:
        The quotient of the ring returned by :func:`order` by
        P^n where P is the prime ideal returned by :func:`prime_ideal`
        and n is the argument.
        
        EXAMPLES::
        
            sage: pAdicBase(ZZ, 3).quotient_ring(2)
            Ring of integers modulo 9
            
        Another example over a number field::
        
            sage: K = QuadraticField(-2)
            sage: pAdics = pAdicBase(K, K.prime_above(7))
            sage: pAdics.quotient_ring(3, names='a')
            Quotient of Maximal Order in Number Field in a with defining polynomial x^2 + 2 by the ideal (343)
            
        """
        if n not in ZZ or n < 0:
            raise ValueError("n should be a non-negative integer, not %s."%(n,))
        if self._R == ZZ:
            return ZZ.quotient((self._P**n).gens()[0])
        else:
            return self._R.quotient(self._P**n, names=names)
    
    def zero(self):
        r"""
        Gives the zero element associated to these p-adics.
        
        OUTPUT:
        The zero element of the order returned by :func:`order`.
        
        EXAMPLE::
        
            sage: pAdicBase(ZZ, 11).zero()
            0
            
        """
        return self._R.zero()
        
    def _repr_(self):
        return "p-adic information of %s with respect to %s."%(
                                                           self.number_field(),
                                                           self._P._repr_short())
                                                           
    def is_extension_of(self, other):
        r"""
        Determines whether these p-adics extend another.
        
        We say this pAdicBase extends another pAdicBase iff
        - The number field of the other pAdicBase is a subfield
          of the number field of this pAdicBase object.
        - The prime ideal of this pAdicBase object lies above
          the prime ideal of the other pAdicBase object.
          
        INPUT:
        
        - ``other`` -- A pAdicBase object.
        
        OUTPUT:
        True if this pAdicBase object extends the pAdicBase object
        given, False in all other cases.
        
        EXAMPLES::
        
            sage: K = QuadraticField(-7)
            sage: p1 = pAdicBase(QQ, 7)
            sage: P = K.prime_above(7)
            sage: p2 = pAdicBase(K, P)
            sage: p2.is_extension_of(p1)
            True
            
        If one field does not extend another::
        
            sage: K = QuadraticField(-3)
            sage: L = QuadraticField(5)
            sage: P = K.prime_above(2)
            sage: Q = L.prime_above(2)
            sage: p1 = pAdicBase(K, P)
            sage: p2 = pAdicBase(L, Q)
            sage: p1.is_extension_of(p2)
            False
        
        If the primes do not match::
        
            sage: K = CyclotomicField(3)
            sage: P = K.prime_above(3)
            sage: p1 = pAdicBase(QQ, 2)
            sage: p2 = pAdicBase(K, P)
            sage: p2.is_extension_of(p1)
            False
        """
        if not isinstance(other, pAdicBase):
            return false
        if not other.number_field().is_subring(self.number_field()):
            return false
        P = other.prime_ideal()
        if other.number_field() == QQ:
            if self.number_field() == QQ:
                return P == self.prime_ideal()
            P = P.gens()[0]
        return self.prime_ideal() in self.number_field().primes_above(P)
        
    def extension_multiplicity(self, other):
        r"""
        Gives the extension multiplicity of some pAdicBase object
        over another.
        
        If some pAdicBase object extends another, then its extension
        multiplicity is the ramification index of the seconds
        defining prime ideal over the ring of integers associated
        to the first.
        
        INPUT:
        
        - ``other`` -- A pAdicBase object that extends this pAdicBase
          object, or of which this pAdicBase object is an extension.
          
        OUTPUT:
        A strictly positive integer equal to the valuation of a
        uniformizer. The valuation will be of the pAdicBase object
        that extends the other and the uniformizer will be from the
        pAdicBase object that is extended.
        
        EXAMPLES::
        
            sage: K = QuadraticField(-5)
            sage: P = K.prime_above(5)
            sage: p1 = pAdicBase(QQ, 5)
            sage: p2 = pAdicBase(K, P)
            sage: p1.extension_multiplicity(p2)
            2
            sage: P.ramification_index()
            2
        
        An example not over QQ::
            
            sage: K = CyclotomicField(4)
            sage: L = CyclotomicField(8)
            sage: P = K.prime_above(3)
            sage: Q = L.prime_above(P)
            sage: p1 = pAdicBase(K, P)
            sage: p2 = pAdicBase(L, Q)
            sage: p1.extension_multiplicity(p2)
            1
        
        It produces the right result when other==self::
        
            sage: K = QuadraticField(3)
            sage: P = K.prime_above(7)
            sage: pAdics = pAdicBase(K, P)
            sage: pAdics.extension_multiplicity(pAdics)
            1
        """
        if self.is_extension_of(other):
            return self.valuation(other.uniformizer())
        elif other.is_extension_of(self):
            return other.valuation(self.uniformizer())
        else:
            raise ValueError("%s and %s do not extend one another."%(self,
                                                                     other))
    
    def _extension_vector_space(self, other):
        K = other.number_field()
        P = other.prime_ideal()
        key = (K,P)
        if key not in self._ext:
            F = self.residue_field()
            L = self.number_field()
            gamma = F.multiplicative_generator()
            VF = F.vector_space()
            G = other.residue_field()
            p = self.characteristic()
            VG = G.vector_space()
            m = VG.rank()
            n = ZZ(VF.rank()/m)
            WG = G^n
            
            e = [G(a) for a in VG.basis()] #Basis of G
            Fe = [F(L(G.lift(ei))) for ei in e] #Embedding of that basis in F
            # Making the matrix that describes the map VG^n -> VF
            # Note the matrix will be on the right of the vector in multiplication!
            MF = Fe
            gammaFe = Fe
            for i in range(1,n):
                gammaFe = [gamma * ei for ei in gammaFe]
                MF.extend(gammaFe)
            M = matrix([VF(MFi._vector_()) for MFi in MF])
            N = M.inverse() #The matrix describing VF -> VG^n
            
            def phi(x):
                x = F(x)
                Vx = VF(x._vector_()) * N
                result = [G(VG(Vx[m*i:m*(i+1)])) for i in range(n)]
                return WG(result)
            
            self._ext[key] = phi
        return self._ext[key]
                                                                     
    def extension_vector_space(self, other):
        """r
        Gives the residue field of this pAdicBase object
        as a vector space over another.
        
        If the other is an extension of this pAdicBase
        object, then it will give the other as a pAdicBase
        over the residue field of this pAdicBase object.
        
        INPUT:
        
        - ``other`` - A pAdicBase object which is an extension
          of this pAdicBase object or which this pAdicBase
          object extends.
          
        OUTPUT:
        A `G`-vectorspace isomorphism `F \to G^n`, where `F` and
        `G` are the residue fields of self and other, with
        `F` being an extension of `G`.
        """
        if self.is_extension_of(other):
            return self._extension_vector_space(other)
        elif other.is_extension_of(self):
            return other.extension_vector_space(self)
        else:
            raise ValueError("%s and %s do not extend one another."%(self,
                                                                     other))
        
    def __eq__(self, other):
        return isinstance(other, pAdicBase) \
               and other.order() == self.order() \
               and other.prime_ideal() == self.prime_ideal()
           
    def __ne__(self, other):
        return not isinstance(other, pAdicBase) \
               or other.order() != self.order() \
               or other.prime_ideal() != self.prime_ideal()
