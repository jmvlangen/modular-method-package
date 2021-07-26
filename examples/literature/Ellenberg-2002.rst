=========================================================================================================
 Galois representations attached to Q-curves and the generalized Fermat equation :math:`A^4 + B^2 = C^p`
=========================================================================================================

We run the computations for some results in the article "Galois
representations attached to Q-curves and the generalized Fermat
equation :math:`A^4 + B^2 = C^p`" written by Jordan S. Ellenberg and
published in American Journal of Mathematics, volume 126 (2004),
no. 4, pages 763-787.

.. linkall

This article considers solutions to the equation :math:`A^4 + B^2 =
C^p` for :math:`A`, :math:`B` and :math:`C` coprime integers and
:math:`p` a prime number at least 211. We will use the notation ``Cp``
to denote :math:`C^p`.

::

   sage: load('load.sage')
   sage: R.<A, B> = QQ[]
   sage: Cp = A^4 + B^2


To each non-trivial solution of the above equation, i.e. :math:`AB \ne
0`, the article attaches a Frey curve that is also a Q-curve. The
article here assumes without loss of generality that :math:`B` is not
congruent to 1 modulo 4.

::

   sage: con = (CoprimeCondition([A, B]) & ~CongruenceCondition(B-1, 4) &
   ....:        PowerCondition(Cp, 3))
   sage: K.<i> = QuadraticField(-1)
   sage: a_invariants = [0, 2*(1+i)*A, 0, B + i*A^2, 0]
   sage: E = FreyQcurve(a_invariants, condition=con, guessed_degrees=[2])

We check that the invariants :math:`c_4` and :math:`\Delta` as in the
article are as mentioned.

::

   sage: E.c4() == 80*i*A^2 - 48*B
   True
   sage: E.discriminant() == -64*i*(A^2 + i*B)*(A^2 - i*B)^2
   True

The article assumes that the restriction of scalars of this curve over
:math:`\QQ(i)` is an abelian surface of GL_2-type. We verify that
:math:`\QQ(i)` is the decomposition field of the curve and that it
decomposes in precisely one factor. We also verify that this single
factor has real multiplication by :math:`\sqrt{2}`.

::

   sage: E.decomposition_field().is_isomorphic(K)
   True
   sage: E.does_decompose()
   True
   sage: E.splitting_image_field('conjugacy')
   (Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?,)

As claimed in proposition 4.5 in the article the level of an
associated newform is indeed 32 or 256. Here we use that the only
possible bad prime is the prime above 2.

::

   sage: E.newform_levels(bad_primes=K.primes_above(2))
   [(256,)] if ('A', 'B') == (1, 0) mod 2
   [(32,)]  if ('A', 'B') == (0, 3), (2, 3) mod 4

The article now uses properties not available in the framework to
prove that these newforms can not be associated with the Frey curve
when :math:`p \ge 211`.
