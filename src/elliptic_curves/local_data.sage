r"""An instance of EllipticCurveLocalData that works for Frey curves.

AUTHORS:

- Joey van Langen (2020-05-15): version ported from frey_curves.sage
  to prevent circular reference

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol
from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class

class FreyCurveLocalData(EllipticCurveLocalData):
    r"""The class for the local reduction data of a Frey curve.

    """
    def __init__(self, E, P, conductor_valuation,
                 discriminant_valuation, kodaira_symbol,
                 tamagawa_number, reduction_type, urst, urst_inv):
        r"""Initialize the reduction data for the Frey curve `E` at the
        prime `P`.

        INPUT:

        - ``E`` -- a Frey curve.

        - ``P`` -- a prime of the definition field of `E`, given as a
          prime number if the definition field is $\QQ$ or as a prime
          ideal otherwise.

        - ``conductor_valuation`` -- The valuation of the conductor of
          `E` at the prime `P`.

        - ``discriminant_valuation`` -- The valuation of the
          discriminant of `E` at the prime `P`.

        - ``kodaira_symbol`` -- The Kodaira symbol associated to the
          reduction of `E` at the prime `P`

        - ``tamagawa_number`` -- The Tamagawa number associated to the
          reduction of `E` at the prime `P`.

        - ``reduction_type`` -- The reduction type of `E` at the prime
          `P`, which can be the values: None, for good reduction, 1
          for split multiplicative reduction, -1 for non-split
          multiplicative reduction, and 0 for additive reduction.

        - ``urst`` -- The change in Weierstrass model to get from the
          original curve to the minimal model stored in this object.

        - ``urst_inv`` -- The change in Weierstrass model to get from
          the minimal model stored in this object to the original
          curve.

        """
        self._set_elliptic_curve(E)
        self._prime=P
        self._fp = conductor_valuation
        self._val_disc = discriminant_valuation
        if isinstance(kodaira_symbol, KodairaSymbol_class):
            self._KS = kodaira_symbol
        else:
            self._KS = KodairaSymbol(kodaira_symbol)
        self._cp = tamagawa_number
        self._reduction_type = reduction_type
        self._urst = urst
        self._urst_inv = urst_inv
        
    def _set_elliptic_curve(self, elliptic_curve):
        self._Emin = elliptic_curve
        self._Emin_reduced = elliptic_curve
        
    def same_local_model(self, other):
        r"""Tell if the reduction data of this object and another are the same.

        INPUT:

        - ``other`` -- Some object

        OUTPUT:

        True, if `other` is an instance of
        :class:`EllipticCurveLocalData` and contains the same
        reduction data as this object. False, otherwise

        """
        return (isinstance(other, EllipticCurveLocalData) and
                self.prime() == other.prime() and
                self.kodaira_symbol() == other.kodaira_symbol() and
                self.conductor_valuation() == other.conductor_valuation() and
                self.tamagawa_number() == other.tamagawa_number())
        
    def same_elliptic_data(self, other):
        r"""Tell if the data of this object and another are the
        same.

        INPUT:

        - ``other`` -- Some object

        OUTPUT:

        True, if `other` is an instance of
        :class:`EllipticCurveLocalData` and contains the same data as
        this object. False, otherwise

        """
        return (self.same_local_model(other) and
                self.minimal_model == other.minimal_model())
    
    def __eq__(self, other):
        return (isinstance(other, ParametrizedLocalData) and
                self.same_elliptic_data(other))
        
    def __ne__(self, other):
        return (not isinstance(other, ParametrizedLocalData) or
                not self.same_elliptic_data(other))
