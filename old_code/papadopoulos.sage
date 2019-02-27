r"""
Unfinished work: An attempt at computing the conductor of a Frey curve
using the table of Papadopolous
"""

def _calculate_vC4(case

def _papadopoulus_invariant_cases(E, T, verbose=False, precision_cap=20):
    S = T.polynomial_ring()
    

def _papadopoulus_lookup_char_p(E, T, verbose=False, precision_cap=20,
                                restrictions=[]):
    
    
    return result

def perform_papadopoulus_lookup((elliptic_curve, initial_values=None,
    prime=None, coefficient_ring=None, base_ring=None, coPrimality=0,
    verbose=False, precision_cap=20, only_calculate=[]):
    r"""
    Calculates local data by doing a table lookup as provided in the
    article by Papadopoulus.
    
    INPUT:
    
    - ``elliptic_curve`` -- An elliptic curve of which the
      coefficients of its Weierstrass equation are inside the
      polynomial ring coefficient_ring. Tate's algorithm will be
      performed on this elliptic curve.
      
    - ``initial_values`` -- A pAdicTree (default: None) containing
      the possible values for the parameters in the given
      elliptic curve. If not given, this will be initialized by
      this function, in which case both the arguments prime and
      coefficient_ring must be defined.
      
    - ``prime`` -- A prime
      (default: initial_values.pAdics().prime_ideal()) inside
      ``base_ring`` for which Tate's algorithm must be performed.
      This may be given as an ideal in base_ring or a generator
      thereof.
      
    - ``coefficient_ring`` -- A polynomial ring over base_ring
      (default: initial_values.polynomial_ring()) of which the
      variables are precisely the parameters of elliptic_curve.
      
    - ``base_ring`` -- (default: coefficient_ring.base_ring()) A
      number field, or subring thereof, on which Tate's algorithm
      can be performed. If not specified it will be the coefficient
      ring of coefficient_ring which only works if coefficient_ring
      is actually a polynomial ring. In most cases this will be the
      ring of integers of a number field.
      
    - ``coPrimality`` -- A non-negative integer (default: 0).
      This option can be set to an integer n if it is known
      that the parameters are n-wise coprime, where n is the
      value of this parameter. The initial values for the
      parameters will then be set up such that no n parameters
      are simultaneously divisible by the given prime. If an
      initial case was set with the argument ``initial_values``
      then this argument is ignored. This argument can at most
      be the number of parameters.
      
    - ``verbose`` -- A boolean value (default: False). If set
      to True the program will print information about the steps
      it is performing and will ask its subfunctions to do the same.
      
    - ``precision_cap`` -- A non-negative integer (default: 20).
      This argument determines the highest precision that will be
      used for the approximation of the parameters. Note that
      setting this too low might result in inaccurate results, for
      which a warning will appear if this is the cases. Setting
      this argument too high might result in endless and slow
      computations.
    
    - ``only_calculate`` -- (default: []) When set to a
      specific value or list of values, will ensure that the
      algorithm only calculates a certain quantity or quantities.
      Possible values are:
        - 'conductor' > Only calculate the conductor exponent
        - 'discriminant' > Only calculate the valuation of the
                           discriminant
        - 'type' > Only calculate the Kodaira Symbol
        - 'minimal_model' > Only calculate the minimal model
                            for this elliptic curve
      Note that some calculations might require some other
      quantities to be determined as well, so the reduction in
      computation might depend on what choices are necessary.
    
    OUTPUT:
    
    The output is a list of possible outcomes of Tate's algorithm
    when performed on the elliptic curve ``elliptic_curve`` for
    the specified prime. When the argument only_calculate is not
    specified, each entry in this list is an object of the type
    ParametrizedLocalData containing the local information about
    this object, combined with the necessary values of the
    parameters stored in a pAdicTree. If however only_calculate
    was specified each entry will be a list wherein each entry
    corresponds to the requested value in only_calculate and
    the last entry is the pAdicTree of value for which these
    values are achieved.
    
    EXAMPLES:
    
    To be made
    
    """
    if initial_values is None:
        if coefficient_ring is None or prime is None:
            raise ValueError("If no initial values are given, both the arguments coefficient_ring and prime must be given.")
        if base_ring is None:
            base_ring = coefficient_ring.base_ring()
        variables = list(coefficient_ring.variable_names())
        pAdics = pAdicBase(base_ring, prime, width=len(variables))
        initial_values = pAdicTree(root=pAdicNode(pAdics=pAdics, full=True),
                                variables=variables)
        if coPrimality not in ZZ or coPrimality < 0 or coPrimality > pAdics.width:
            raise ValueError("The argument coPrimality is not well-defined")
        if coPrimality > 0:
            initial_values.apply_coprimality_restriction(coPrimality=coPrimality)
    if isinstance(only_calculate, str):
        only_calculate = [only_calculate]
    
    pAdics = initial_values.pAdics()
    char = pAdics.characteristic()
    if char == 2:
        l = pAdics.valuation(2)
        if l == 1:
            return _papadopoulus_lookup_char2_val1(elliptic_curve,
                                                   initial_values,
                                                   verbose=verbose,
                                                   precision_cap=precision_cap,
                                                   restrictions=only_calculate)
        else:
            return _papadopoulus_lookup_char2(elliptic_curve,
                                              initial_values, l
                                              verbose=verbose,
                                              precision_cap=precision_cap,
                                              restrictions=only_calculate)
    elif char == 3:
        l = pAdics.valuation(3)
        if l == 1:
            return _papadopoulus_lookup_char3_val1(elliptic_curve,
                                                   initial_values,
                                                   verbose=verbose,
                                                   precision_cap=precision_cap,
                                                   restrictions=only_calculate)
        else:
            return _papadopoulus_lookup_char3(elliptic_curve,
                                              initial_values, l
                                              verbose=verbose,
                                              precision_cap=precision_cap,
                                              restrictions=only_calculate)     
    else:
        return _papadopoulus_lookup_char_p(elliptic_curve,
                                           initial_values,
                                           verbose=verbose,
                                           precision_cap=precision_cap,
                                           restrictions=only_calculate)
        
