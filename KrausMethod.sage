def containsTwoZeros(iterableObject):
    count = 0
    for item in iterableObject:
        if item == 0:
            count += 1
        if count >= 2:
            return True
    return False

# Creates a list of all possible a_p values for an elliptic curve defined with parameters over QQ.
# ellipticCurve : An elliptic curve, dependant on some parameters.
# primeRange : A list of prime numbers for which possible a_p values should be listed.
#       This could also be an integer in which case the list will be all primes up to that integer.
# [primesToAvoid] : A list of prime numbers that should be excluded from the range.
#       By default this list is empty.
# [printProgress] : If set to True the program will print progress mesaages during execution.
#       By default this is set to False.
# [pairwiseCoprimality] : If set to True the program will exclude possibilities for which the variables are definitely not pairwise coprime..
#       By default this is set to False.
# return value : A list with for each prime in primeRange that is not in primesToAvoid a list containing as first entry the corresponding prime p and as second entry a list of all corresponding values a_p that can be optained by substituting values for the parameters in the elliptic curve ellipticCurve.
def getApList( ellipticCurve , primeRange , primesToAvoid = [] , printProgress = False , pairwiseCoprimality = False):
    ApList = []
    
    # Initializing the primes
    if type(primeRange) is not list:
        primeRange = prime_range(primeRange)
    for p in primesToAvoid:
        if p in primeRange:
            primeRange.remove(p)
    
    # Retrieving the relevant variables
    variables , numberOfVariables = getVariables( list(ellipticCurve.a_invariants()) )
    
    # Retrieving the discriminant
    delta = symbolic_expression(ellipticCurve.discriminant())
    delta = delta.expand()

    if printProgress :
        print "Initialization complete: parameters" , variables

    # The algorithm
    ApList = []
    for p in primeRange:
        Fp = FiniteField(p)
        values = Fp^numberOfVariables
        S = []
        for value in values:
            if not (pairwiseCoprimality and containsTwoZeros(value)):
                evaluatedDelta = delta
                for i in range(numberOfVariables):
                    evaluatedDelta = evaluatedDelta.subs(variables[i] == value[i])
                evaluatedDelta = Fp(evaluatedDelta)
                if evaluatedDelta != 0:
                    evaluatedAinvariants = copy(aInvariants)
                    for j in range(len(aInvariants)):
                        for i in range(numberOfVariables):
                            evaluatedAinvariants[j] = evaluatedAinvariants[j].subs(variables[i] == value[i])
                        evaluatedAinvariants[j] = Fp(evaluatedAinvariants[j])
                    Ee = EllipticCurve(evaluatedAinvariants)
                    ap = p + 1 - Ee.count_points()
                    if ap not in S:
                        S.append(ap)
        ApList.append([p,S])
        if printProgress:
            print "Calculated case: p =" , p
        
    return ApList

# Applies the method of Kraus to see which newforms can match a given list of possible values a_p.
# This method basicly looks whether the mod l Galois representation of an elliptic curve and a newform can match.
# ApList : A list consististing of lists with as first entry a prime number p and as second entry a list containing possible values of a_p for some elliptic curve.
# newForms : A list of newforms that should be compared to the given ApList
# [printProgress] : Can be set to True if the code should print progress messages during execution.
# result value : A list consisting of lists with as first entry a newform f which could match the list of given values a_p and as second entry the possible primes l for which the Galois representations could match. If it works for all l then this second entry does not exist.   
def applyKrausMethod( ApList , newForms , printProgress = False):
    result = []
    
    for f in newForms:
        N = f.level()
        B = 0
        for [p,S] in ApList:
            if not p.divides(N):
                af = f.coefficient(p)
                if af.base_ring == QQ:
                    Bp = p * ((p + 1)^2 - af^2)
                else:
                    Bp = p * norm( (p+1)^2 - af^2 )
                for ap in S:
                    if Bp == 0:
                        break
                    af = f.coefficient(p)
                    if af.base_ring == QQ:
                        Bp = Bp * af - ap
                    else:
                        Bp = Bp * norm(af - ap)
                B = ZZ(gcd(B,Bp))
                if B == 1:
                    break
        if printProgress:
            print "Calculated case: f =" , f , "B =" , B
        if B == 0:
            result.append([f])
        else:
            if B != 1:
                result.append([f,B.prime_factors()])
      
    return result
    
