def getVariables( listOfExpressions ):
    variables = []
    for i in range(len(listOfExpressions)):
        listOfExpressions[i] = symbolic_expression(listOfExpressions[i])
        if sage.symbolic.ring.is_SymbolicExpressionRing(listOfExpressions[i].base_ring()):
            for v in listOfExpressions[i].variables():
                if not v in variables:
                    variables.append(v)
    return variables, len(variables)

def getEllipticCurveFromBinaryCubicForm( binaryCubicForm ):
    F = binaryCubicForm
    if not F.is_homogeneous():
        raise ValueError("Form is not homogeneous")
    if F.degree() != 3:
        raise ValueError("Form is not cubic")
    
    # Detecting the variables
    variables = list(F.variables())
    if len(variables) != 2:
        raise ValueError("Form is not binary")
    x = variables[0]
    y = variables[1]
    
    # Calculating necessary invariants
    H = -1 / 4 * ( F.derivative(x).derivative(x) * F.derivative(y).derivative(y) - (F.derivative(x).derivative(y))^2 )
    G = F.derivative(x) * H.derivative(y) - F.derivative(y) * H.derivative(x)
    Delta = F(x,1).discriminant(x)
    
    # Check whether everything is alright
    if not 4 * H^3 == G^2 +  27 * Delta * F^2:
        raise ValueError("Something went wrong whilst checking some constants")
    
    return EllipticCurve([-3*H,G])
   
def indicesNotSimultaneouslyZero( m , S ):
    for s in S:
        flag = True
        for i in s:
            if m[i] != 0:
                flag = False
                break
        if flag:
            return False
    return True
    
 
def getCorrespondingZquotient( ring ):
    f = ZZ.hom(ring)
    n = 1
    while f(n) != 0:
        n += 1
    return ZZ.quotient(n)
    
def getUniformizer( primeIdeal , ring ):
    if ring == ZZ:
        return primeIdeal.gens()[0]
    for g in primeIdeal.gens():
        if g.ord( primeIdeal ) == 1:
            return g
    raise ValueError("Somehow this ideal does not have a uniformizer")

def getAllExtensions( element , Zmodule0 , Zmodule ):
    N0 = Zmodule0.base_ring().cardinality()
    N = Zmodule.base_ring().cardinality()
    r = Zmodule0.rank()
    k = N / N0
    indexModule = ( ZZ.quotient(k) )^r
    result = []
    for m in indexModule:
        result.append( Zmodule( element ) + Zmodule( m ) )
    return result

def getIdealChain( RidealStart , RidealEnd , R ):
    if R == ZZ:
        n0 = RidealStart.gens()[0]
        n1 = RidealEnd.gens()[0]
        k = (n1/n0).factor()
        totalFactors = sum( [ f[1] for f in k ] )
        
        result = [ ZZ.ideal[n0] ]
        for j in range(totalFactors):
            j0 = j
            for f in k:
                if j < f[1]:
                    result.append(result[-1] * ZZ.ideal(f[0]))
                    break
                else:
                    j0 -= f[1]
        return result
    else:
        k = (RidealEnd / RidealStart).factor()
        totalFactors = sum( [ f[1] for f in k ] )
        result = [ RidealStart ]
        for j in range(totalFactors):
            j0 = j
            for f in k:
                if j < f[1]:
                    result.append(result[-1] * f[0])
                    break
                else:
                    j0 -= f[1]
        return result

def getModuleChains( RidealStart , RidealEnd , R , numberOfVariables ):
    idealChain = getIdealChain( RidealStart , RidealEnd , R )
    RmoduleChain = []
    ZmoduleChain = []
    for I in idealChain:
        RmodI = R.quotient( I )
        ZmodI = getCorrespondingZquotient( RmodI )
        RmoduleChain.append( RmodI^numberOfVariables )
        ZmoduleChain.append( ZmodI^numberOfVariables )
    
    return ZmoduleChain , RmoduleChain
    

def getZeroCasesFromModuleChains( element , index , ZmoduleChain , RmoduleChain , invariant ):
    Zmodule0 = ZmoduleChain[index - 1]
    Zmodule = ZmoduleChain[index]
    Rmodule = RmoduleChain[index]
    
    zeroCases = []
    nonZeroCases = []
    for e in getAllExtensions( element , Zmodule0 , Zmodule ):
        
        if index < len( ZmoduleChain ) - 1:
            if invariant( tuple( Rmodule( e ) ) ) == 0:
                newZeroCases , newNonZeroCase = getZeroCasesFromModuleChains( e , index + 1 , rank , ZmoduleChain , RmoduleChain , invariant )
                zeroCases.extend( newZeroCases )
                nonZeroCases.extend( newNonZeroCases )
            else:
                nonZeroCases.extend( getAllExtensions( e , Zmodule , ZmoduleChain[-1] ) )
        
        else:
            if invariant( tuple( Rmodule( e ) ) ) == 0:
                zeroCases.append( e )
            else:
                nonZeroCases.append( e )
    
    return zeroCases , nonZeroCases
    
def getZeroCases( M0 , vals , R , P , k , invariant ):
    RmoduleChain = []
    ZmoduleChain = []
    r = M0.rank()
    for i in range(k+1):
        Rmod = R.quotient( P , names = 'b' )
        Zmod = getCorrespondingZquotient( Rmod )
        RmoduleChain.append( Rmod^r )
        ZmoduleChain.append( Zmod^r )
    while ZmoduleChain[0] != M0:
        RmoduleChain.pop(0)
        ZmoduleChain.pop(0)
           
    zeroCases = []
    nonZeroCase = []
    if len(RmoduleChain) == 1:
        for e in vals:
            if invariant( tuple( RmoduleChain[0]( e ) ) ) == 0:
                zeroCases.append(e)
            else:
                nonZeroCases.append(e)
    else:
        for e in vals:
            newZeroCases , newNonZeroCases = getZeroCasesFromModuleChains( e , 1 , ZmoduleChain , RmoduleChain , invariant )
            zeroCases.extend(newZeroCases)
            nonZeroCases.extend(newNonZeroCases)
    
    return zeroCases , nonZeroCases , ZmoduleChain[-1]
    
def getCasesInWhichInvariantFitsInPr( E , P , R , S , variables, numberOfVariables, M0 , vals, invariant,  r , caseYesFirstEntry , caseNoFirstEntry ):
    if invariant in R:
        if invariant in P^r:
            return [ [ caseYesFirstEntry , [ M0 , vals ] , E ] ]
        else:
            return [ [ caseNoFirstEntry , [ M0 , vals ] , E ] ]
    else:
        # Determining power of P already present without substituting parameters.
        if R == ZZ:
            n = r - gcd(invariant.coefficients()).ord(P.gens()[0])
        else:
            n = r - gcd(invariant.coefficients()).ord(P)
        if n <= 0:
            return [ [ caseYesFirstEntry , [ M0 , vals ] , E ] ]
            
        RmodP = R.quotient(P^r,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M2 = ZmodP^numberOfVariables
        N = RmodP^numberOfVariables
        RmodP = R.quotient(P^n,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M1 = ZmodP^numberOfVariables

        # Determining best space for parameters.
        if M1.cardinality() > M0.cardinality():
            M = M1
        else:
            M = M0

        caseYes = [ caseYesFirstEntry , [ M , [] ] , E ]
        caseNo = [ caseNoFirstEntry , [ M , [] ] , E ]
        result = []
        
        for m in M:
            if M0( m ) in vals:
                if invariant( tuple( N( M2( m ) ) ) ) == 0:
                    caseYes[1][1].append(m)
                else:
                    caseNo[1][1].append(m)

        if len( caseYes[1][1] ) > 0:
            result.append(caseYes)
        if len( caseNo[1][1] ) > 0:
            result.append(caseNo)
                    
        return result

def getMaximalOrderAmongCoefficients( invariant , prime , ring ):
    if invariant == 0:
        return Infinity
    if ring == ZZ:
        return gcd(invariant.coefficients()).ord( prime.gens()[0] )
    return gcd( invariant.coefficients() ).ord( prime )
    
    
def performTatesAlgorithmStep1( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( E.discriminant() ),
        1,
        2,
        [ "Type I0" , 0 , 1 , 0 , 1 ] )
    
def performTatesAlgorithmStep2Transformation( E , P , R , S , variables , numberOfVariables , M0 , vals , char ):
    RmodP = R.quotient(P,names='b')
    ZmodP = getCorrespondingZquotient(RmodP)
    M1 = ZmodP^numberOfVariables
    N = RmodP^numberOfVariables
    N2 = RmodP^2
    replaceCases = dict()
    result = []
    
    # Determining best space for parameters.
    if M1.cardinality() > M0.cardinality():
        M = M1
    else:
        M = M0
    
    # First transforming the singular point to (0,0), which (could) give us different cases.
    s = 0
    if char == 2:
        if R == ZZ and P in ZZ.ideal_monoid():
            s = ( ( ZZ.quotient( P.gen() - 1 )(2) )^(-1) ).lift()
        else:
            s = ( ( ZZ.quotient( ZZ(norm(P) - 1) )(2) )^(-1) ).lift()
    if char == 3:
        if R == ZZ and P in ZZ.ideal_monoid():
            s = ( ( ZZ.quotient( P.gen() - 1 )(3) )^(-1) ).lift() 
        else:
            s = ( ( ZZ.quotient( ZZ( norm(P) - 1 ) )(3) )^(-1) ).lift()
    if s == 0:
        s = 1
    
    for m in M:
        if (M0(m) in vals):
            # Isolating the singular point x,y
            a1 = S(E.a1())( tuple( N( M1( m ) ) ) )
            a2 = S(E.a2())( tuple( N( M1( m ) ) ) )
            a3 = S(E.a3())( tuple( N( M1( m ) ) ) )
            a4 = S(E.a4())( tuple( N( M1( m ) ) ) )
            a6 = S(E.a6())( tuple( N( M1( m ) ) ) )
        
            if char == 2:
                if a1 == 0:
                    x = a4^s
                    y = (a6 + a2 * a4)^s
                else:
                    x = a4 / a1
                    y = (x^2 + a4) / a1
            else:
                # Coordinate transformation to end up in the case where a1 = a3 = 0,
                # hence y = 0
                a22 = a2 + a1^2 / 4
                a42 = a4 + a1 * a3 / 2
                a62 = a6 + a3^2 / 4
                y = 0
                
                if char == 3:
                    if a22 == 0:
                        x = (-a62)^s
                    else:
                        x = - a42 / a22
                else:
                    # Coordinate transformation to end up in the case where also a2 = 0
                    a43 = a42 - a22^2 / 3
                    a63 = a62 - a22 * a42 / 3 + 2 * a22^3 / 27
                    
                    if a43 == 0:
                        x = 0
                    else:
                        x = - 3 * a63 / ( 2 * a43 )
                    
                    # Transforming back
                    x = x - a22 / 3
                    
                # Transforming back
                y = y - a1 * x / 2 - a3 / 2
            
            singularPoint = tuple( N2( [x,y] ) )
            if singularPoint in replaceCases:
                replaceCases[singularPoint].append(m)
            else:
                replaceCases[singularPoint] = [m]
    
    for (point,values) in replaceCases.iteritems():
        x0 = point[0].lift()
        y0 = point[1].lift()
        E0 = E.rst_transform(x0,0,y0)
        result.append( [ 2 + 1/2 , [ M , values ] , E0 ] )
            
    return result
    
def performTatesAlgorithmStep2( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( E.b2() ),
        1,
        3,
        [ "Type In" , "v( Delta )" , "v( Delta )" , 1 , [ 1 , 2 , "v( Delta )" ] ] )

def performTatesAlgorithmStep3( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( E.a6() ),
        2,
        4,
        [ "Type II" , "v( Delta )" , 1 , "v( Delta )" , 1 ] )
    
def performTatesAlgorithmStep4( E , P , R , S , variables , numberOfVariables , M0 , vals ):   
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( E.b8() ),
        3,
        5,
        [ "Type III" , "v( Delta )" , 2 , "v( Delta ) - 1" , 2 ] )

def performTatesAlgorithmStep5( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( E.b6() ),
        3,
        6,
        [ "Type IV" , "v( Delta )" , 3 , "v( Delta ) - 2" , [3,1] ] )

def performTatesAlgorithmStep6Transformation( E , P , R , S , variables , numberOfVariables , M0 , vals , char ):
    pi = getUniformizer( P , R )
    if char != 2:
    
        a1 = S( E.a1() )
        a3 = S( E.a3() )
        
        k = max( 1 - getMaximalOrderAmongCoefficients( a1 , P , R ) ,
            2 - getMaximalOrderAmongCoefficients( a3 , P , R ) )
        if k <= 0:
            return [ [ 6 + 1/2 , [ M0 , vals ] , E ] ]
            
        RmodP = R.quotient( P , names='b' )
        ZmodP = getCorrespondingZquotient(RmodP)
        N1 = RmodP^numberOfVariables
        M1 = ZmodP^numberOfVariables      
        half = ( RmodP( 2 )^(-1) )
        
        # Determining best space for parameters.
        if M1.cardinality() > M0.cardinality():
            M = M1
        else:
            M = M0
        
        changeDict = dict()    
        for m in M:
            if M0( m ) in vals:
                alpha = - half * a1( tuple( N1( M( m ) ) ) )
                a31 = RmodP( R( a3( tuple( ( N1( M( m ) ) ).lift() ) ) / pi ) )
                beta = -half * a31
                alphaBetaPair = tuple( [ alpha, beta ] )
                if alphaBetaPair in changeDict:
                    changeDict[alphaBetaPair].append( m )
                else:
                    changeDict[alphaBetaPair] = [ m ]
        
        result = []
        for ( alphaBetaPair , values ) in changeDict.iteritems():
            E0 = E.rst_transform( 0 , ( M.base_ring()( alphaBetaPair[0] ) ).lift() , ( M.base_ring()( alphaBetaPair[1] ) ).lift() * pi  )
            result.append( [ 6 + 1/2 , [ M , values ] , E0 ] )
        return result
    else:
        if R == ZZ and P in ZZ.ideal_monoid():
            sqrtPower = ( ( ZZ.quotient( P.gen() - 1 )(2) )^(-1) ).lift()
        else:
            sqrtPower = ( ( ZZ.quotient( ZZ(norm(P) - 1) )(2) )^(-1) ).lift()
        if sqrtPower == 0:
            sqrtPower = 1
            
        a2 = S( E.a2() )
        a6 = S( E.a6() )
        k1 = 3 - getMaximalOrderAmongCoefficients( a6 , P , R )
        k2 = 1 - getMaximalOrderAmongCoefficients( a2 , P , R )
        k = max(k1,k2)
        if k <= 0:
            return [ [ 6 + 1/2 , [ M0 , vals ] , E ] ]
        
        RmodP = R.quotient(P^3,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        Ma6 = ZmodP^numberOfVariables
        Na6 = RmodP^numberOfVariables
        RmodP = R.quotient(P^k,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M1 = ZmodP^numberOfVariables
        RmodP = R.quotient(P , names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        Ma2 = ZmodP^numberOfVariables 
        Na2 = RmodP^numberOfVariables

        # Determining best space for parameters.
        if M1.cardinality() > M0.cardinality():
            M = M1
        else:
            M = M0
        
        changeDict = dict()    
        for m in M:
            if M0( m ) in vals:
                alphaSqrd = - a2( tuple( Na2 ( Ma2 ( m ) ) ) )
                betaSqrd = RmodP( ( - a6( tuple( Na6 ( Ma6 ( m ) ) ) ) ).lift() / pi^2 )
                alpha = alphaSqrd^sqrtPower
                beta = betaSqrd^sqrtPower
                alphaBetaPair = tuple( [ alpha, beta ] )
                if alphaBetaPair in changeDict:
                    changeDict[alphaBetaPair].append( m )
                else:
                    changeDict[alphaBetaPair] = [ m ]
        
        result = []
        for ( alphaBetaPair , values ) in changeDict.iteritems():
            E0 = E.rst_transform( 0 , ( M.base_ring()( alphaBetaPair[0] ) ).lift() , ( M.base_ring()( alphaBetaPair[1] ) ).lift() * pi  )
            result.append( [ 6 + 1/2 , [ M , values ] , E0 ] )
        return result

def performTatesAlgorithmStep6( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( -4 * E.a2()^3 * E.a6() + E.a2()^2 * E.a4()^2 - 4 * E.a4()^3 - 27 * E.a6()^2 + 18 * E.a2() * E.a4() * E.a6() ),
        7,
        7,
        [ "Type I0*" , "v( Delta )" , 5 , "v( Delta ) - 4" , [2,4] ] )
        
def performTatesAlgorithmStep7( E , P , R , S , variables , numberOfVariables , M0 , vals , char ):
    if char == 2:
        return getCasesInWhichInvariantFitsInPr(
            E, P, R, S, variables, numberOfVariables, M0, vals,
            S( 3 * E.a4() - ( E.a2() )^2 ),
            3,
            8,
            7 + encodeQuotient( 1 , 0 ) )
    else:
        return getCasesInWhichInvariantFitsInPr(
            E, P, R, S, variables, numberOfVariables, M0, vals,
            S( 3 * E.a4() - ( E.a2() )^2 ),
            3,
            8,
            [ "Type In*" , "v( Delta )" , "v( Delta ) + 5" , "v( Delta ) - 4" , [2,4] ] )

def encodeQuotient( n , b ):
    m = 2 * n + b
    return 1 / ( m + 1 )

def decodeQuotient( q ):
    m = ZZ( 1 / q ) - 1
    if is_even(m):
        n = ZZ( m / 2 )
        b = 0
    else:
        n = ZZ( ( m - 1 ) / 2 )
        b = 1   
    return n , b
            
def performTatesAlgorithmStep7subTransformation( E , P , R , S , variables , numberOfVariables , M0 , vals , n ):
    if R == ZZ and P in ZZ.ideal_monoid():
        sqrtPower = ( ( ZZ.quotient( P.gen() - 1 )(2) )^(-1) ).lift()
    else:
        sqrtPower = ( ( ZZ.quotient( ZZ(norm(P) - 1) )(2) )^(-1) ).lift()
    if sqrtPower == 0:
        sqrtPower = 1
        
    if n == 1:
        RmodP = R.quotient( P )
        k2 = 3 - getMaximalOrderAmongCoefficients( S( E.a4() ) , P , R )
        if k2 <= 0:
            return [ [ 7 + encodeQuotient( n , 1 ) , [ M0 , vals ] , E ] ]
            
        RmodP = R.quotient(P^3,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M3 = ZmodP^numberOfVariables
        N3 = RmodP^numberOfVariables
        RmodP = R.quotient(P^k2,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M1 = ZmodP^numberOfVariables
        pi = getUniformizer(P , R)
        RmodP = R.quotient(P , names='b')

        # Determining best space for parameters.
        if M1.cardinality() > M0.cardinality():
            M = M1
        else:
            M = M0
        
        changeDict = dict()    
        for m in M:
            if M0( m ) in vals:
                a4 = ( S( E.a4() )( tuple( N3 ( M3 ( m ) ) ) ) ).lift()
                square = RmodP( R( a4 / ( pi^2 ) ) )
                change = square^sqrtPower
                if change in changeDict:
                    changeDict[change].append( m )
                else:
                    changeDict[change] = [ m ]
        
        result = []
        for ( change , values ) in changeDict.iteritems():
            E0 = E.rst_transform( pi * change.lift() , 0 , 0 )
            result.append( [ 7 + encodeQuotient( n , 1 ) , [ M , values ] , E0 ] )
        
        return result
    else:
        k = 3 + n
        if R == ZZ:
            if E.a6() == 0:
                k2 = -Infinity
            else:
                k2 = k - gcd(S( E.a6() ).coefficients()).ord(P.gens()[0])
            if is_odd( n ) and E.a2() != 0:
                k2 = max( k2 , 2 - gcd(S( E.a2() ).coefficients()).ord(P.gens()[0]) )
        else:
            if E.a6() == 0:
                k2 = -Infinity
            else:
                k2 = k - gcd(S( E.a6() ).coefficients()).ord(P)
            if is_odd( n ) and E.a2() != 0:
                k2 = max( k2 , 2 - gcd(S( E.a2() ).coefficients()).ord(P.gens()[0]) )
        if k2 <= 0:
            return [ [ 7 + encodeQuotient( n , 1 ) , [ M0 , vals ] , E ] ]
        
        RmodP = R.quotient(P^k,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        Mk = ZmodP^numberOfVariables
        Nk = RmodP^numberOfVariables
        RmodP = R.quotient(P^2,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M2 = ZmodP^numberOfVariables
        N2 = RmodP^numberOfVariables
        RmodP = R.quotient(P^k2,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M1 = ZmodP^numberOfVariables
        pi = getUniformizer(P , R)
        RmodP = R.quotient(P , names='b')

        # Determining best space for parameters.
        if M1.cardinality() > M0.cardinality():
            M = M1
        else:
            M = M0
        
        changeDict = dict()    
        for m in M:
            if M0( m ) in vals:
                a6 = ( S( E.a6() )( tuple( Nk ( Mk ( m ) ) ) ) ).lift()
                a2 = ( S( E.a2() )( tuple( N2 ( M2 ( m ) ) ) ) ).lift()
                square = RmodP( R( a6 / ( pi^(k-1) ) ) )
                if is_odd(n):
                    square = square / ( RmodP( R( a2 / pi ) ) )
                change = square^sqrtPower
                if change in changeDict:
                    changeDict[change].append( m )
                else:
                    changeDict[change] = [ m ]
        
        result = []
        for ( change , values ) in changeDict.iteritems():
            if is_odd( n ):
                E0 = E.rst_transform( - pi^( ZZ( ( n + 1 ) / 2 ) ) * change.lift() , 0 , 0 )
            else:
                E0 = E.rst_transform( 0 , 0 , pi^( ZZ( ( n + 2 ) / 2 ) ) * change.lift() )
            result.append( [ 7 + encodeQuotient( n , 1 ) , [ M , values ] , E0 ] )
        
        return result
        
def performTatesAlgorithmStep7sub( E , P , R , S , variables , numberOfVariables , M0 , vals , n ):
    if is_odd(n):
        return getCasesInWhichInvariantFitsInPr(
            E, P, R, S, variables, numberOfVariables, M0, vals,
            S( E.a3() ),
            ZZ( 3 + (n-1)/2 ),
            7 + encodeQuotient( n + 1 , 0 ),
            [ "Type I" + str( n ) + "*" , "v( Delta )" , 5 + n , "v( Delta ) - " + str( 4 + n ) , [2,4] ] )
    else:
        return getCasesInWhichInvariantFitsInPr(
            E, P, R, S, variables, numberOfVariables, M0, vals,
            S( E.a4() ),
            ZZ( 3 + n/2 ),
            7 + encodeQuotient( n + 1 , 0 ),
            [ "Type I" + str( n ) + "*" , "v( Delta )" , 5 + n , "v( Delta ) - " + str( 4 + n ) , [2,4] ] )
     
def performTatesAlgorithmStep8Transformation( E , P , R , S , variables , numberOfVariables , M0 , vals , char ):
    if char != 3:
        RmodP = R.quotient( P , names='b' )
        return [ [ 8 + 1/2 , [ M0 , vals ] , E.rst_transform( - E.a2() * ( RmodP( 3 )^(-1) ).lift() , 0 , 0 ) ] ]
    else:
        if R == ZZ and P in ZZ.ideal_monoid():
            cubertPower = ( ( ZZ.quotient( P.gen() - 1 )(3) )^(-1) ).lift()
        else:
            cubertPower = ( ( ZZ.quotient( ZZ(norm(P) - 1) )(3) )^(-1) ).lift()
        k = 4 - getMaximalOrderAmongCoefficients( S( E.a6() ) , P , R )
        if k <= 0:
            return [ [ 8 + 1/2 , [ M0 , vals ] , E ] ]
        
        RmodP = R.quotient(P^4,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M4 = ZmodP^numberOfVariables
        N4 = RmodP^numberOfVariables
        RmodP = R.quotient(P^k,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M1 = ZmodP^numberOfVariables
        pi = getUniformizer(P , R)
        RmodP = R.quotient(P , names='b')

        # Determining best space for parameters.
        if M1.cardinality() > M0.cardinality():
            M = M1
        else:
            M = M0
        
        changeDict = dict()    
        for m in M:
            if M0( m ) in vals:
                a6 = ( S( E.a6() )( tuple( N4 ( M4 ( m ) ) ) ) ).lift()
                cube = RmodP( R( a6 / ( pi^3 ) ) )
                change = cube^cubertPower
                if change in changeDict:
                    changeDict[change].append( m )
                else:
                    changeDict[change] = [ m ]
        
        result = []
        for ( change , values ) in changeDict.iteritems():
            E0 = E.rst_transform( -pi * change.lift() , 0, 0  )
            result.append( [ 8 + 1/2 , [ M , values ] , E0 ] )
        
        return result
            
def performTatesAlgorithmStep8( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( ( E.a3() )^2 + 4 * E.a6() ),
        5,
        9,
        [ "Type IV*" , "v( Delta )" , 7 , "v( Delta ) - 6" , [3,1] ] )
        
def performTatesAlgorithmStep9Transformation( E , P , R , S , variables , numberOfVariables , M0 , vals , char ):
    if char != 2:
        RmodP = R.quotient( P , names='b' )
        return [ [ 9 + 1/2 , [ M0 , vals ] , E.rst_transform( 0 , 0 , - E.a3() * ( RmodP( 2 )^(-1) ).lift() ) ] ]
    else:
        if R == ZZ and P in ZZ.ideal_monoid():
            sqrtPower = ( ( ZZ.quotient( P.gen() - 1 )(3) )^(-1) ).lift()
        else:
            sqrtPower = ( ( ZZ.quotient( ZZ(norm(P) - 1) )(3) )^(-1) ).lift()
        if sqrtPower == 0:
            sqrtPower = 1
        if E.a6() == 0:
            k = -Infinity
        elif R == ZZ:
            k = 5 - gcd(S( E.a6() ).coefficients()).ord(P.gens()[0])
        else:
            k = 5 - gcd(S( E.a6() ).coefficients()).ord(P)
        if k <= 0:
            return [ [ 9 + 1/2 , [ M0 , vals ] , E ] ]
        
        RmodP = R.quotient(P^5,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M5 = ZmodP^numberOfVariables
        N5 = RmodP^numberOfVariables
        RmodP = R.quotient(P^k,names='b')
        ZmodP = getCorrespondingZquotient(RmodP)
        M1 = ZmodP^numberOfVariables
        pi = getUniformizer(P , R)
        RmodP = R.quotient(P , names='b')

        # Determining best space for parameters.
        if M1.cardinality() > M0.cardinality():
            M = M1
        else:
            M = M0
        
        changeDict = dict()    
        for m in M:
            if M0( m ) in vals:
                a6 = ( S( E.a6() )( tuple( N5 ( M5 ( m ) ) ) ) ).lift()
                square = RmodP( R( -a6 / ( pi^4 ) ) )
                change = square^sqrtPower
                if change in changeDict:
                    changeDict[change].append( m )
                else:
                    changeDict[change] = [ m ]
        
        result = []
        for ( change , values ) in changeDict.iteritems():
            E0 = E.rst_transform( 0, 0, - pi^2 * change.lift()  )
            result.append( [ 9 + 1/2 , [ M , values ] , E0 ] )
        
        return result
            
def performTatesAlgorithmStep9( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( E.a4() ),
        4,
        10,
        [ "Type III*" , "v( Delta )" , 8 , "v( Delta ) - 7" , 2 ] )
            
def performTatesAlgorithmStep10( E , P , R , S , variables , numberOfVariables , M0 , vals ):
    return getCasesInWhichInvariantFitsInPr(
        E, P, R, S, variables, numberOfVariables, M0, vals,
        S( E.a6() ),
        6,
        "Not minimal",
        [ "Type II*" , "v( Delta )" , 9 , "v( Delta ) - 8" , 1 ] )
 
def moduloCaseCleanUp( moduloCase ):
    M = moduloCase[0]
    vals = moduloCase[1]
    r = M.rank()
    Zmod = M.base_ring()
    N = Zmod.cardinality()
    
    for d in N.divisors():
        if (d^r).divides( len( vals ) ) and d != 1:
            M1 = ( ZZ.quotient( ZZ( N / d ) ) )^r
            vals1 = []
            for e in vals:
                if M1 ( e ) not in vals1:
                    vals1.append( M1( e ) )
            if len( vals1 ) * d^r == len( vals ):
                return moduloCaseCleanUp( [ M1 , vals1 ] )
    
    return [ M , vals ]
    
def casesCleanUp( cases ):
    for case in cases:
        case[1] = moduloCaseCleanUp( case[1] )
    
def performTatesAlgorithm( ellipticCurve , prime , coefficientRing ,
    baseRing=0 , initialCase = 0 , coPrimality = 0,  printProgress=False):
    r"""
    Performs Tate's Algorithm on an elliptic curve which may depend on certain parameters.
    
    In case the curve depends on parameters, it will make a case distinction
    with for each case the relevant parameter values (modulo a suitable prime) the 
    elliptic curve will fit into that case.
    For now all possible parameter values are considered to be integers.
    
    INPUT:
    
    - ``ellipticCurve`` -- An elliptic curve of which the coefficients of its Weierstrass equation
      are inside the ring coefficientRing. Tate's algorithm will be performed on this elliptic curve.
    - ``prime`` -- A prime inside ``baseRing`` for which Tate's algorithm must be performed.
      This may be given as an ideal in baseRing or a generator thereof.
    - ``coefficientRing`` -- A polynomial ring over baseRing of which the variables are precisely
      the parameters of ellipticCurve. If no parameters must be considered then coefficientRing
      must equal baseRing.
    - ``baseRing`` -- (default: coefficientRing.base_ring()) A Dedekind domain on which Tate's algorithm can be performed.
      If not specified it will be the coefficient ring of coefficientRing which only works
      if coefficientRing is actually a polynomial ring in at least one variable.
      In most cases this will be the ring of integers of a number field.
    - ``initialCase`` -- (default: None) The case in which the program should execute.
      This should be formatted as entry [1] of a case as described under Output.
      i.e. it should be a list with as first entry a free module over a quotient of ZZ
      with the standard basis elements corresponding to the respective parameters,
      and the second entry is a list of values in this module which the parameters could attain.
    - ``coPrimality`` -- (default: 0) This option can be set to an integer n if it is known
      that the parameters are n-wise coprime, where n is the value of this parameter.
      The initial case will then be set up such that no n parameters are simultaneously
      divisible by the prime ``prime``. If an initial case was set with the argument ``initialCase``
      then this argument is ignored.
      This argument can at most be the number of parameters and must at least be 0.
    - ``printProgress`` -- (default: False) If set to True the program will print the step it is performing at the moment.
    
    OUTPUT:
    
    The output is a list of possible cases that could be the outcome of Tate's algorithm
    when performed on the elliptic curve ``ellipticCurve`` for the specified prime.
    Each case is itself a list with as entries:
    
    0. A list detailing the reduction information that results from Tate's algorithm in this
       case, including
       
       0. The reduction type of the curve as expressed in Kodaira symbols.
       1. The number of components over an algebraic closure of the reduction field.
       2. The valuation of the minimal discriminant.
       3. The exponent of the conductor.
       4. The order of the group of components. During execution this list may also just be
          a number, in which case this number is the next step of Tate's algorithm to be
          performed on this case.
         
    1. A list containing information about for which parameter values we are in this case.
    
       0. A free module over a quotient of ZZ of which the standard basis elements
          correspond to the parameters of ellipticCurve.
       1. A list of values which the parameters could take inside the free module for
          which we end up in this case. Note that a tuple of parameter values (x0 , ... , xn )
          makes the elliptic curve end up in this case if (x0 , ... , xn) seen as an element of
          the module 0. ends up in the list 1.
          
    2. The elliptic curve in this case. This might be different from the original elliptic
       curve as a change of variables might be needed for minimality or to ease computation.
    
    EXAMPLES:
    
    To be made
    
    """
    if baseRing == 0:
        baseRing = coefficientRing.base_ring()
    # Easy renaming for programming, must replace back afterwards!
    R = baseRing
    S = coefficientRing
    if prime in R:
        P = R.ideal(prime)
    else:
        P = prime
        if P.number_field() != R.number_field():
            raise ValueError("The prime is not in the coefficientring")
    char = getCorrespondingZquotient( R.quotient(P,names='b') ).characteristic()
    variables = list(S.gens())
    if printProgress:
        print "Detected variables: " , variables
    numberOfVariables = len(variables)
    if initialCase == 0:
        if coPrimality not in ZZ or coPrimality < 0 or coPrimality > numberOfVariables:
            raise ValueError("The argument coPrimality is not well-defined")
        if coPrimality > 0:
            coPrimeCollection = Set( range( numberOfVariables ) ).subsets(coPrimality)
            ZmodChar = ZZ.quotient(char)
            M0 = ZmodChar^numberOfVariables
            newCases = [ [ 1 , [ M0 , [ m for m in M0 if indicesNotSimultaneouslyZero(m,coPrimeCollection) ] ] , ellipticCurve ] ] 
        else:
            M0 = (ZZ.quotient(1))^numberOfVariables
            newCases = [ [ 1 , [ M0 , [M0.zero()] ] , ellipticCurve ] ]
    else:
        newCases = [ [ 1 , initialCase , ellipticCurve ] ]
    doneCases = []
    
    while len(newCases) > 0:
        cases = newCases
        newCases = []
        for case in cases:
            if case[0] == 1:
                if printProgress:
                    print "Performing Step 1"
                newCases.extend( performTatesAlgorithmStep1( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 2:         
                if printProgress:
                    print "Performing Step 2 transformation"  
                newCases.extend( performTatesAlgorithmStep2Transformation( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] , char ) )
            elif case[0] == 2 + 1/2:
                if printProgress:
                    print "Performing Step 2"
                newCases.extend( performTatesAlgorithmStep2( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 3:
                if printProgress:
                    print "Performing Step 3"
                newCases.extend( performTatesAlgorithmStep3( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 4:
                if printProgress:
                    print "Performing Step 4"
                newCases.extend( performTatesAlgorithmStep4( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 5:
                if printProgress:
                    print "Performing Step 5"
                newCases.extend( performTatesAlgorithmStep5( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 6:
                if printProgress:
                    print "Performing Step 6 transformation"
                newCases.extend( performTatesAlgorithmStep6Transformation( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] , char ) )
            elif case[0] == 6 + 1/2:
                if printProgress:
                    print "Performing Step 6"
                newCases.extend( performTatesAlgorithmStep6( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 7:
                if printProgress:
                    print "Performing Step 7"
                newCases.extend( performTatesAlgorithmStep7( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] , char ) )
            elif case[0] in QQ and case[0] < 8:
                n , b = decodeQuotient( case[0] - 7 )
                if b == 0:
                    if printProgress:
                        print "Performing Step 7sub transformation"
                    newCases.extend( performTatesAlgorithmStep7subTransformation( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] , n ) )
                else:
                    if printProgress:
                        print "Performing Step 7sub"
                    newCases.extend( performTatesAlgorithmStep7sub( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] , n ) )
            elif case[0] == 8:
                if printProgress:
                    print "Performing Step 8 transformation"
                newCases.extend( performTatesAlgorithmStep8Transformation( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] , char ) )
            elif case[0] == 8 + 1/2:
                if printProgress:
                    print "Performing Step 8"
                newCases.extend( performTatesAlgorithmStep8( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 9:
                if printProgress:
                    print "Performing Step 9 transformation"
                newCases.extend( performTatesAlgorithmStep9Transformation( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] , char ) )
            elif case[0] == 9 + 1/2:
                if printProgress:
                    print "Performing Step 9"
                newCases.extend( performTatesAlgorithmStep9( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            elif case[0] == 10:
                if printProgress:
                    print "Performing Step 10"
                newCases.extend( performTatesAlgorithmStep10( case[2] , P , R , S , variables , numberOfVariables , case[1][0] , case[1][1] ) )
            else:
                doneCases.append(case)
        
        casesCleanUp( newCases )
        
    casesCleanUp( doneCases )
    return doneCases
    
def getPossibleCandidatesModulo( polynomial , n , variablesOfNoInterest=[]):
    variables, numberOfVariables = getVariables( [ polynomial ] )
    Zmodn = ZZ.quotient(n)
    valueCandidates = Zmodn^numberOfVariables
    
    variablesOfInterest = []
    for v in variables:
        if v not in variablesOfNoInterest:
            variablesOfInterest.append(v)
    relevantValueSpace = Zmodn^len( variablesOfInterest )
    
    relevantValues = []
    for value in valueCandidates:
        if polynomial( tuple( value ) ) == 0:
            relevantValue = []
            for i in range( numberOfVariables ):
                if variables[i] in variablesOfInterest:
                    relevantValue.append( value[i] )
            relevantValues.append( relevantValueSpace( relevantValue ) )
    
    return [ relevantValueSpace , relevantValues ]
