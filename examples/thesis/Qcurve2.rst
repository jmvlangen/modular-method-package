========================
 Code for Example 2.3.1
========================

This file contains the code for Example 2.3.1.

.. linkall

The following import is required for the example to work

::

   sage: from modular_method import *

We verify the computations given in the example.

::

   sage: K.<sqrt2> = QuadraticField(2)
   sage: E = Qcurve([0, 12*sqrt2, 0, 36*(1 + sqrt2), 0], guessed_degrees=[2])
   sage: E.complete_definition_field() == K
   True
   sage: G = K.galois_group()
   sage: [E.degree_map(s) for s in G]
   [1, 2]
   sage: matrix([[E.c(s, t) for t in G] for s in G])
   [ 1  1]
   [ 1 -2]

Next we show that the map given is indeed a splitting map

::

   sage: L.<sqrtm2> = QuadraticField(-2)
   sage: beta = {G[0]: 1, G[1]: sqrtm2}
   sage: all(E.c(s, t) == beta[s]*beta[t]*beta[s*t]^(-1)
   ....:     for s in G for t in G)
   True
