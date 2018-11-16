import itertools

def eliminate_newforms_by_trace(curves, newforms, condition=None, primes=50,
                                precision_cap=1, verbose=False):
    r"""
    Eliminates newforms associated to frey curves.

    Given one or multiple frey curves with the same
    parameters, compares traces of frobenius at different
    primes to see which newforms can and which newforms
    can not be associated to these curves.

    INPUT:

    - ``curves`` -- A frey curve or a list/tuple of frey
      curves that share the same parameters.
    - ``newforms`` -- A list of newforms that are
      associated to the given frey curve, i.e. such that
      their mod-l representations could be isomorphic
      to the mod-l representation of the frey curve.
      Instead of a list one can also give the value
      None, in which case the list will be retrieved
      from the curve by the method
      :meth:`newform_candidates` of the given Frey curve.
      If multiple Frey curves are given this argument
      should be a tuple of length equal to the number
      of given Frey curves, containing for each Frey
      curve an entry as described above in the same order
      as the Frey curves are given.
      Alternatively the input may also be a list of
      tuples wherein each tuple has length equal to
      the number of Frey curves and each entry thereof
      is a newform associated to the corresponding Frey
      curve.
      In the last input format, each tuple may also be
      1 longer than the number of Frey curves, in which
      case the last entry must be an integer divisible
      by all primes l for which the mod-l representations
      of the given elliptic curves may still be
      isomorphic to the corresponding newforms in this
      tuple. This can be used to chain multiple
      elimination methods.
      Everywhere where a list is expected, one may also
      give a ConditionalValue of which the possible
      values are lists of the expected format. The
      method will be applied to each entry individually
      and the condition will be replaced by the condition
      that both the given condition and the condition
      at which that value is attained hold.
    - ``condition`` -- A Condition giving the
      restrictions on the parameters of the given
      Frey curve that should be considered. By default
      this will be set to the condition stored in
      the first given Frey curve.
    - ``primes`` -- A list of prime numbers or a
      strictly positive integer (default: 50). This
      list gives all the primes at which the traces
      of frobenius of the different galois
      representations should be compared. If set to
      strictly positive integer, will be initialized
      as the list of all prime numbers less than
      the given number.
    - ``precision_cap`` -- A non-negative integer
      (default: 1) giving the maximal precision used
      to compute the reduction type of a Frey curve
      at a prime. Note that setting this value to
      a value higher than 1 might result in a long
      computation time and heavy memory usage and
      is strongly discouraged.
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
    
    A list of tuples, wherein each tuple contains:
    - for each frey curve given a corresponding newform
      for which an isomorphism between the mod-l
      representations of that curve and that newform
      could not be excluded by this method.
    - as its last entry an integer that is divisible
      by all prime numbers l for which the mentioned
      mod-l galois representations could still exist.
    """
    if isinstance(curves, FreyCurve):
        curves = (curves,)
    else:
        curves = tuple(curve for curve in curves)
        parameters = None
        if len(curves) == 0:
            raise ValueError("At least one curve should be provided.")            
        for curve in curves:
            if not isinstance(curve, FreyCurve):
                raise ValueError("%s is not a Frey curve"%(curve,))
            if parameters is None:
                parameters = curve.parameters()
            elif parameters != curve.parameters():
                raise ValueError("%s and %s don't have the same parameters"%(curves[0], curve))
    newforms = _init_newform_list(newforms, curves)
    if primes in ZZ and primes > 0:
        primes = prime_range(primes)
    if condition is None:
        condition = curves[0]._condition
    for prime in primes:
        if verbose > 0:
            print "Comparing traces of frobenius at %s for %s cases."%(prime, len(newforms))
        newforms = apply_to_conditional_value(lambda nfs, con:
                                              _eliminate_newforms_by_trace(curves,
                                                                           nfs,
                                                                           con & condition,
                                                                           prime,
                                                                           precision_cap,
                                                                           verbose),
                                              newforms,
                                              use_condition=True,
                                              default_condition=condition)
    return newforms
        
def _eliminate_newforms_by_trace(curves, newforms, condition, prime,
                                 precision_cap, verbose):
    r"""
    Eliminates newforms associated to frey curves.

    Does what is required for a single prime.
    
    see :meth:`eliminate_newforms_by_trace`
    """
    if len(newforms) == 0:
        return newforms
    nE = len(curves)
    fields = tuple(curve.definition_field() for curve in curves)
    primes = tuple((prime if K == QQ else K.prime_above(prime)) for K in fields)
    powers = tuple((1 if P in ZZ else P.residue_class_degree() * P.ramification_index()) for P in primes)
    apE_ls = _init_traces(curves, condition, primes, precision_cap, verbose)
    result = []
    for nfs in newforms:
        Bold = nfs[-1]
        if all(not prime.divides(nfs[i].level()) for i in range(nE)):
            apf = [nfs[i].trace_of_frobenius(prime, power=powers[i]) for i in range(nE)]
            B = ZZ(prime * product(gcd([(QQ(apE[i] - apf[i]) if apf[i] in QQ
                                         else (apE[i] - apf[i]).absolute_norm())
                                        for i in range(nE)])
                                   for apE in apE_ls))
            Bnew = gcd(Bold, B)
        else:
            Bnew = Bold
        if abs(Bnew) != 1:
            result.append(tuple([nfs[i] for i in range(nE)] + [Bnew]))
    return result

def _init_traces(curves, condition, primes, precision_cap, verbose):
    eP = [(1 if prime in ZZ else prime.ramification_index()) for prime in primes]
    traces = [curves[i].trace_of_frobenius(primes[i], power=eP[i],
                                           condition=condition,
                                           precision_cap=precision_cap,
                                           verbose=(verbose - 1 if verbose > 0 else verbose))
              for i in range(len(curves))]
    pAdics = pAdicBase(curves[0].definition_ring(), primes[0]).pAdics_below(curves[0]._R)
    default_tree = condition.pAdic_tree(pAdics=pAdics,
                                        verbose=(max(0, verbose - 3) if verbose > 0 else verbose),
                                        precision_cap=precision_cap)
    traces = [(trace if isinstance(trace, ConditionalValue)
               else [(trace, TreeCondition(default_tree))])
              for trace in traces]
    result = []
    for case in itertools.product(*traces):
        values, conditions = zip(*case)
        trees = [con.pAdic_tree() for con in conditions]
        tree = trees[0]
        for i in range(1, len(trees)):
            tree = tree.intersection(trees[i])
        if not tree.is_empty():
            result.append(values)
    return result
                                       
def _init_newform_list(newforms, curves):
    """
    Initializes a list of newforms associated
    to n frey curves
    """
    if isinstance(newforms, tuple):
        if len(newforms) != len(curves):
            raise ValueError("Expected %s lists of newforms, but got %s"%(len(curves),
                                                                          len(newforms)))
        newforms = [(apply_to_conditional_value(lambda x: list(x),
                                                curves[i].newform_candidates())
                     if newforms[i] is None
                     else newforms[i])
                    for i in range(len(curves))]
        newforms = conditional_product(*newforms)
        newforms = apply_to_conditional_value(lambda nfs: itertools.product(*nfs),
                                              newforms)
    if isinstance(newforms, ConditionalValue):
        return apply_to_conditional_value(lambda nfs: _init_newform_list(nfs, curves), newforms)
    newforms = [nfs for nfs in newforms]
    for i in range(len(newforms)):
        if not isinstance(newforms[i], tuple):
            if len(curves) != 1:
                raise ValueError("Expected %s newforms per entry, but got 1"%(len(curves),))
            newforms[i] = (newforms[i], ZZ(0))
        if len(newforms[i]) == len(curves) + 1:
            pass
        elif len(newforms[i]) == len(curves):
            newforms[i] = tuple(list(newforms[i]) + [0])
        else:
            raise ValueError("Expected %s newforms per entry, but got %s"%(len(curves),
                                                                           len(newforms[i])))
    return newforms

def kraus_method(curves, newforms, l, polynomials, primes=200,
                 condition=None, precision_cap=1, verbose=False):
    r"""
    Eliminates newforms associated to Frey curves with
    the Kraus method.

    Given one or multiple frey curves with the same
    parameters, a prime number l and a list of polynomials
    in the variables that are l-th powers, computes
    which of the given newforms can still be associated
    when comparing traces of frobenius of primes which
    are 1 mod l for their representations mod l.

    INPUT:

    - ``curves`` -- A frey curve or a list/tuple of frey
      curves that share the same parameters.
    - ``newforms`` -- A list of newforms that are
      associated to the given frey curve, i.e. such that
      their mod-l representations could be isomorphic
      to the mod-l representation of the frey curve.
      Instead of a list one can also give the value
      None, in which case the list will be retrieved
      from the curve by the method
      :meth:`newform_candidates` of the given Frey curve.
      If multiple Frey curves are given this argument
      should be a tuple of length equal to the number
      of given Frey curves, containing for each Frey
      curve an entry as described above in the same order
      as the Frey curves are given.
      Alternatively the input may also be a list of
      tuples wherein each tuple has length equal to
      the number of Frey curves and each entry thereof
      is a newform associated to the corresponding Frey
      curve.
      In the last input format, each tuple may also be
      1 longer than the number of Frey curves, in which
      case the last entry must be an integer divisible
      by all primes l for which the mod-l representations
      of the given elliptic curves may still be
      isomorphic to the corresponding newforms in this
      tuple. This can be used to chain multiple
      elimination methods.
      Everywhere where a list is expected, one may also
      give a ConditionalValue of which the possible
      values are lists of the expected format. The
      method will be applied to each entry individually
      and the condition will be replaced by the condition
      that both the given condition and the condition
      at which that value is attained hold.
    - ``prime`` -- A prime number.
    - ``polynomials`` -- A polynomial or list thereof
      in the common parameters of the given Frey curves,
      such that these polynomials are l-th powers.
    - ``primes`` -- A list of prime numbers or a
      strictly positive integer (default: 200). This
      list should contain all prime numbers for which
      the traces of frobenius at that prime should be
      compared. If it is set to a strictly positive
      integer will initialize this as the list
      prime_range(l, primes). Note that only those
      primes for which the relevant residue fields
      have order congruent to 1 mod l will be
      actually considered.
    - ``condition`` -- A Condition giving the
      restrictions on the parameters of the given
      Frey curve that should be considered. By default
      this will be set to the condition stored in
      the first given Frey curve.
    - ``precision_cap`` -- A non-negative integer
      (default: 1) giving the maximal precision used
      to compute the reduction type of a Frey curve
      at a prime. Note that setting this value to
      a value higher than 1 might result in a long
      computation time and heavy memory usage and
      is strongly discouraged.
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
    
    A list of tuples, wherein each tuple contains:
    - for each frey curve given a corresponding newform
      for which an isomorphism between the mod-l
      representations of that curve and that newform
      could not be excluded by this method.
    - as its last entry an integer that is divisible
      by all prime numbers l for which the mentioned
      mod-l galois representations could still exist.
    """
    if isinstance(curves, FreyCurve):
        curves = (curves,)
    else:
        curves = tuple(curve for curve in curves)
        parameters = None
        if len(curves) == 0:
            raise ValueError("At least one curve should be provided.")            
        for curve in curves:
            if not isinstance(curve, FreyCurve):
                raise ValueError("%s is not a Frey curve"%(curve,))
            if parameters is None:
                parameters = curve.parameters()
            elif parameters != curve.parameters():
                raise ValueError("%s and %s don't have the same parameters"%(curves[0], curve))
    newforms = _init_newform_list(newforms, curves)
    if not (l in ZZ and l > 0 and l.is_prime()):
        raise ValueError("%s is not a prime number"%(prime,))
    if not isinstance(polynomials, list):
        if hasattr(polynomials, "__iter__"):
            polynomials = [poly for poly in polynomials]
        else:
            polynomials = [polynomials]
    if primes in ZZ and primes > 0:
        primes = prime_range(l, primes)
    if not (isinstance(primes, list) and all(p in ZZ and p > 0 and p.is_prime() for p in primes)):
        raise ValueError("%s is not a list of prime numbers"%(primes,))
    if condition is None:
        condition = curves[0]._condition
        
    fields = tuple(poly.base_ring() for poly in polynomials)
    for p in primes:
        Ps = [_primes_1_mod_l(K, p, l) for K in fields]
        for P in itertools.product(*Ps):
            CP = [_init_kraus_condition(polynomials[i], P[i], l) for i in range(len(P))]
            if verbose > 0:
                print "Comparing traces of frobenius at %s for %s cases."%(p, len(newforms))
            C = reduce(lambda x, y: x & y, C, condition)
            newforms = apply_to_conditional_value(lambda nfs, con:
                                                  _eliminate_newforms_by_trace(curves,
                                                                               nfs,
                                                                               l,
                                                                               prime,
                                                                               con & C,
                                                                               precision_cap,
                                                                               verbose),
                                                  newforms,
                                                  use_condition=True,
                                                  default_condition=C)
    return newforms

def _primes_1_mod_l(K, p, l):
    if K == QQ:
        if mod(p, l) == 1:
            return [p]
        else:
            return []
    else:
        return [P for P in K.primes_above(p)
                if mod(P.residue_field().cardinality(), l) == 1]
                        
def _init_kraus_condition(polynomial, prime, l):
    if prime in ZZ:
        F = Integers(prime)
        u = F.unit_gens()[0]
    else:
        F = prime.residue_field()
        u = F.primitive_element()
    values = [F.lift(0)] + [F.lift(u^(n*l)) for n in range((F.cardinality() - 1)/l)]
    conditions = [CongruenceCondition(polynomial - val, prime) for val in values]
    def combine_conditions(C1, C2):
        if C1 is None:
            return C2
        elif C2 is None:
            return C1
        else:
            return C1 | C2
    return reduce(combine_conditions, conditions, None)

def _kraus_method_trace_compare(curves, newforms, l, p, condition, precision_cap, verbose):
    if len(newforms) == 0:
        return newforms
    nE = len(curves)
    fields = tuple(curve.definition_field() for curve in curves)
    primes = tuple((p if K == QQ else K.prime_above(p)) for K in fields)
    powers = tuple((1 if P in ZZ else P.residue_class_degree() * P.ramification_index()) for P in primes)
    apE_ls = _init_traces(curves, condition, primes, precision_cap, verbose)
    result = []
    for nfs in newforms:
        Bold = nfs[-1]
        if l.divides(Bold) and all(not p.divides(nfs[i].level()) for i in range(nE)):
            apf = [nfs[i].trace_of_frobenius(p, power=powers[i]) for i in range(nE)]
            B = max(min((QQ(apE[i] - apf[i]) if apf[i] in QQ
                         else (apE[i] - apf[i]).absolute_norm()).ord(l)
                        for i in range(nE))
                    for apE in apE_ls)
            if B > 0:
                Bnew = Bold / (l^(Bold.ord(l)))
            else:
                Bnew = Bold
        else:
            Bnew = Bold
        if abs(Bnew) != 1:
            result.append(tuple([nfs[i] for i in range(nE)] + [Bnew]))
    return result
