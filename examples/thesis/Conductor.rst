=========================
 Examples from Chapter 1
=========================

This document contains the code for the examples from Chapter 1 of the
thesis.

.. linkall

The following import is required for all examples to work

::

   from modular_method import *

Example 1.2.1
-------------

We enter the elliptic curve mentioned in this example.

::

   sage: R.<A, B> = QQ[]
   sage: E = EllipticCurve([0, B - A, 0, -A*B, 0])

We check the discriminant is indeed as mentioned in Step 1.

::

   sage: E.discriminant() == 2^4 * A^2 * B^2 * (A + B)^2
   True

We check the equations from Step 2.

::

   sage: F2 = GF(2)
   sage: S.<x0, y0> = R.change_ring(F2)[]
   sage: f = E.defining_polynomial().change_ring(R.change_ring(F2))(x0, y0, 1)
   sage: f
   x0^3 + (A + B)*x0^2 + y0^2 + A*B*x0
   sage: fx = f.derivative(x0); fx
   x0^2 + A*B
   sage: fy = f.derivative(y0); fy
   0

Next we find all the x0 and y0 that satisfy this for all A and B
modulo 2, which agrees with the example.

::

   sage: {(A, B) : [(x0, y0) for x0 in range(2) for y0 in range(2)
   ....:            if (f(x0, y0)(A, B) == 0 and
   ....:                fx(x0, y0)(A, B) == 0 and
   ....:                fy(x0, y0)(A, B) == 0)]
   ....:  for A in range(2) for B in range(2)}
   {(0, 0): [(0, 0)], (0, 1): [(0, 0)], (1, 0): [(0, 0)], (1, 1): [(1, 0)]}

We compute the different models of E that follow from this.

::

   sage: E1 = E.rst_transform(0, 0, 0); E1
   Elliptic Curve defined by y^2 = x^3 + (-A+B)*x^2 + (-A*B)*x over Multivariate Polynomial Ring in A, B over Rational Field
   sage: E2 = E.rst_transform(1, 0, 0); E2
   Elliptic Curve defined by y^2 = x^3 + (-A+B+3)*x^2 + (-A*B-2*A+2*B+3)*x + (-A*B-A+B+1) over Multivariate Polynomial Ring in A, B over Rational Field

We verify that b2 is indeed as mentioned in the example.

::

   sage: E1.b2() == 4*(B - A)
   True
   sage: E2.b2() == 4*(B - A + 3)
   True

Next we check a6 is indeed as in the example.

::

   sage: E1.a6() == 0
   True
   sage: E2.a6() == (1 - A)*(1 + B)
   True

Now we check b8.

::

   sage: E1.b8() == -A^2*B^2
   True
   sage: E2.b8() == 3 - 4*A + 4*B - 6*A*B - A^2*B^2
   True

As described we now compute all possible A and B modulo 8 such that b8
is not divisible by 8, and check this is equal to all A and B such
that v(A B (A + B)) = 1.

::

   sage: b8 = lambda A, B: E1.b8()(A, B) if 2.divides(A*B) else E2.b8()(A, B)
   sage: vb8lt3 = [(A, B) for A in range(8) for B in range(8)
   ....:           if not 8.divides(b8(A, B))]
   sage: vABApBeq1 = [(A, B) for A in range(8) for B in range(8)
   ....:              if (2.divides(A*B*(A + B)) and
   ....:                  not 4.divides(A*B*(A + B)))]
   sage: vb8lt3 == vABApBeq1
   True

Example 1.5.5
-------------

We verify the results from Example 1.5.5 with some computations. First
of all we compute the transform of the curve mentioned.

::

   sage: E7 = E.rst_transform(0, 1, 0); E7
   Elliptic Curve defined by y^2 + 2*x*y = x^3 + (-A+B-1)*x^2 + (-A*B)*x over Multivariate Polynomial Ring in A, B over Rational Field

We now go through the steps of Tate's algorithm. First in step 1 the
discriminant is 2^4 * A^2 * B^2 * (A + B)^2

::

   sage: E7.discriminant() == 2^4 * A^2 * B^2 * (A + B)^2
   True

so we can continue to Step 2. All coefficients of the Weierstrass
equation are divisible by 2, so we can choose the trivial
transformation. Now we compute b2

::

   sage: E7.b2()
   -4*A + 4*B

which is divisible by 2 so we continue to Step 3. As we compute a6 ==
0 we can continue to Step 4 and compute b8.

::

   sage: E7.b8()
   -A^2*B^2

By assumption this is divisible by 2^6 so we continue to Step 5. Now
we compute b6

::

   sage: E7.b6()
   0

so again we can continue. Note that by our earlier conclusions the
right hand side of the equations for \alpha and \beta are zero modulo
2, so we can choose trivial transformations for Step 6. We check that
the invariant mentioned in Step 6 is the discriminant of P(T) times
\pi^6

::

   sage: Ra.<a1, a2, a3, a4, a6, pi> = QQ[]
   sage: Sa.<T> = Ra.fraction_field()[]
   sage: P = T^3 + (a2 / pi)*T^2 + (a4 / pi^2)*T + (a6 / pi^3)
   sage: pi^6 * P.discriminant() == -4*a2^3*a6 + a2^2*a4^2 - 4*a4^3 - 27*a6^2 + 18*a2*a4*a6
   True

To check that the invariant of Step 7 is indeed whether P(T) has a
triple root or not, we start with the assumption that P(T) has a
double root and show this invariant computes whether this root is a
triple root as well.

::

   sage: Rb.<alpha, beta> = QQ[]
   sage: Sb.<T> = Rb[]
   sage: Pb = (T - alpha)^2 * (T - beta)
   sage: a6, a4, a2, a0 = Pb.coefficients()
   sage: (alpha - beta)^2 == a2^2 - 3*a4
   True

The rest is already verified in the example itself.

Example 1.6.1
-------------

The code for this example is the same as given in the example. It is
simply provided here as a way to verify the output.

Firs of all we construct the corresponding FreyCurve object.

::

   sage: from modular_method import *
   sage: R.<A, B> = QQ[]
   sage: con = CoprimeCondition([A, B])
   sage: E = FreyCurve([0, B - A, 0, -A*B, 0], condition=con); E
   Frey curve defined by y^2 = x^3 + (-A+B)*x^2 + (-A*B)*x over Rational Field with parameters (A, B)

Now we can compute the conductor over :math:`\QQ` with one method.
   
::
   
   sage: N = E.conductor(); N
   Warning: Assuming that A and B are coprime.
   2^n0*Rad_P( (16) * B^2 * A^2 * (A + B)^2 )
    where
   n0 = 5 if ('A', 'B') is 1 of 6 possibilities mod 4
        4 if ('A', 'B') is 1 of 3 possibilities mod 4
        3 if ('A', 'B') is 1 of 36 possibilities mod 16
        0 if ('A', 'B') is 1 of 24 possibilities mod 32
        1 if ('A', 'B') is 1 of 24 possibilities mod 32

The warning comes from the necessary assumption that :math:`c_4` and
:math:`\Delta` are coprime outside some finite set of primes. In this
case the finite set was chosen as :math:`{ 2 }` by the default method.
        
::
        
   sage: E.primes_of_possible_additive_reduction()
   [2]

Note that the `Rad_P` part is not explicitly computed. It just
displays the factorisation of the discriminant :math:`\Delta`. It
indicates that the remaining part of the conductor is just the product
of all primes dividing :math:`\Delta` that are not in :math:`\\{ 2
\\}`.

We can also change the set :math:`\\{ 2 \\}` to compute more conductor
exponents explicitly.
   
::
   
   sage: E.conductor(additive_primes=[2, 3, 5, 7])
   2^n0*3^n1*5^n2*7^n3*Rad_P( (16) * B^2 * A^2 * (A + B)^2 )
    where
   n0 =  5 if ('A', 'B') is 1 of 6 possibilities mod 4
         4 if ('A', 'B') is 1 of 3 possibilities mod 4
         3 if ('A', 'B') is 1 of 36 possibilities mod 16
         0 if ('A', 'B') is 1 of 24 possibilities mod 32
         1 if ('A', 'B') is 1 of 24 possibilities mod 32
   n1 =  0 if ('A', 'B') == (1, 1), (2, 2) mod 3
         1 if ('A', 'B') is 1 of 6 possibilities mod 3
   n2 =  0 if ('A', 'B') is 1 of 12 possibilities mod 5
         1 if ('A', 'B') is 1 of 12 possibilities mod 5
   n3 =  0 if ('A', 'B') is 1 of 30 possibilities mod 7
         1 if ('A', 'B') is 1 of 18 possibilities mod 7

We can also impose additional conditions. For example we could compute
the conductor exponent at :math:`2` as in Example 1.5.5.

::

   sage: con2 = (CongruenceCondition(A*B, 8) &
   ....:         CongruenceCondition(A - B - 1, 4))
   sage: E.conductor_exponent(2, condition=con2)
   4

Note that `E` is also the Frey curve associated with Fermat's Last
Theorem in case we take :math:`A`, :math:`B` and :math:`A + B` to be
l-th powers with :math:`l \ge 3`. We also compute the conductor
imposing these additional conditions.

::
   
   sage: conFLT = (con &
   ....:           PowerCondition(A, 3) &
   ....:           PowerCondition(B, 3) &
   ....:           PowerCondition(A + B, 3))
   sage: E.conductor(condition=conFLT)
   2^n0*Rad_P( (16) * B^2 * A^2 * (A + B)^2 )
    where 
   n0 = 3 if ('A', 'B') is 1 of 12 possibilities mod 16
        4 if ('A', 'B') is 1 of 6 possibilities mod 8
        0 if ('A', 'B') is 1 of 24 possibilities mod 32
        1 if ('A', 'B') is 1 of 24 possibilities mod 32

The results of these functions can be conditional values or
conditional expressions. We illustrate how one can inspect such values
with the conductor computed earlier
   
::

   sage: N.left()
   2^n0
    where
   n0 = 5 if ('A', 'B') is 1 of 6 possibilities mod 4
        4 if ('A', 'B') is 1 of 3 possibilities mod 4
        3 if ('A', 'B') is 1 of 36 possibilities mod 16
        0 if ('A', 'B') is 1 of 24 possibilities mod 32
        1 if ('A', 'B') is 1 of 24 possibilities mod 32
   sage: N.left().right()
   5 if ('A', 'B') is 1 of 6 possibilities mod 4
   4 if ('A', 'B') is 1 of 3 possibilities mod 4
   3 if ('A', 'B') is 1 of 36 possibilities mod 16
   0 if ('A', 'B') is 1 of 24 possibilities mod 32
   1 if ('A', 'B') is 1 of 24 possibilities mod 32
   sage: N.left().right()[1]
   (4, The condition that ('A', 'B') == (0, 3), (1, 0), (3, 1) mod 4)
