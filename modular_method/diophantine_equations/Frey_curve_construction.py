def get_symmetric_terms_ring(polynomial):
    R = polynomial.parent();
    variables = R.variable_names()
    if len(variables) != Integer(2) :
        raise ValueError("Can only work with 2 variables.")
    sumSquares = variables[Integer(0) ] + '2' + variables[Integer(1) ] + '2'
    prod = variables[Integer(0) ] + variables[Integer(1) ]
    S = PolynomialRing(R.base(),[sumSquares,prod])
    return S

def polynomial_in_symmetric_terms(polynomial):
    R = get_symmetric_terms_ring(polynomial)
    a, b = polynomial.variables()
    sumSquares, prod = R.gens()
    K = R.base()
    f = polynomial
    result = R(Integer(0) )
    while f != Integer(0) :
        e_ls = f.exponents()
        i = Integer(0) 
        for j in range(len(e_ls)):
            if e_ls[j][Integer(0) ] > e_ls[i][Integer(0) ] \
            or (e_ls[j][Integer(0) ] == e_ls[i][Integer(0) ] and e_ls[j][Integer(1) ] > e_ls[i][Integer(1) ]):
                 i = j
        c = K(f.coefficient(list(e_ls[i])))
        e2 = e_ls[i][Integer(1) ]
        e1 = ZZ((e_ls[i][Integer(0) ] - e2)/Integer(2) )
        result = result + c * sumSquares**e1 * prod**e2
        f = f - c * (a**Integer(2)  + b**Integer(2) )**e1 * (a*b)**e2
    return result
    
def getEllipticCurveFromBinaryCubicForm( binaryCubicForm ):
    F = binaryCubicForm
    if not F.is_homogeneous():
        raise ValueError("Form is not homogeneous")
    if F.degree() != Integer(3) :
        raise ValueError("Form is not cubic")
    
    # Detecting the variables
    variables = list(F.variables())
    if len(variables) != Integer(2) :
        raise ValueError("Form is not binary")
    x = variables[Integer(0) ]
    y = variables[Integer(1) ]
    R = F.base_ring()
    
    # Calculating necessary invariants
    H = -Integer(1)  / Integer(4)  * ( F.derivative(x).derivative(x) * F.derivative(y).derivative(y) - (F.derivative(x).derivative(y))**Integer(2)  )
    G = F.derivative(x) * H.derivative(y) - F.derivative(y) * H.derivative(x)
    Delta = R[x](F(x,Integer(1) )).discriminant()
    
    # Check whether everything is alright
    if not Integer(4)  * H**Integer(3)  == G**Integer(2)  +  Integer(27)  * Delta * F**Integer(2) :
        raise ValueError("Something went wrong whilst checking some constants")
    
    return EllipticCurve([-Integer(3) *H,G])
    
def getEllipticCurveFromBinarySymmetricSixForm(sixForm):
    F = sixForm
    if not F.is_homogeneous():
        raise ValueError("Form is not homogeneous")
    if F.degree() != Integer(6) :
        raise ValueError("Form is not of degree 6")
        
    variables = list(F.variables())
    if len(variables) != Integer(2) :
        raise ValueError("Form is not binary")
    x = variables[Integer(0) ]
    y = variables[Integer(1) ]
    S = F.base_ring()[x,y]
    F0= S(F)
    if F0 != F0(y,x) or F0 != F0(-x,-y):
        raise ValueError("Form is not symmetric") 
    
    G = polynomial_in_symmetric_terms(F)
    E = getEllipticCurveFromBinaryCubicForm(G)
    cfs = E.a_invariants()
    cfs = [c(x**Integer(2)  + y**Integer(2) , x*y) for c in cfs]
    return EllipticCurve(cfs)
    
def my_factor(x):
    if x == Integer(0) :
        print(Integer(0) )
        return
    f = x.factor()
    if f.unit() == Integer(1) :
        result = ""
    else:
        result = str(f.unit().factor())
    for (a, n) in f:
        result += " * (" + str(a) + ")"
        if n != Integer(1) :
            result += "^" + str(n)
    print(result)

def make_degree_dict(poly_list):
    result = {}
    for f in poly_list:
        d = f.degree()
        if d in result:
            result[d].append(f)
        else:
            result[d] = [f]
    return result

def find_degree_n_combinations(poly_dict, n, min_deg=Integer(1) ):
    if n < min_deg:
        if n == Integer(0) :
            return [Integer(1) ]
        else:
            return []
    else:
        result = []
        for d in range(min_deg, n+Integer(1) ):
            if d in poly_dict:
                g_list = find_degree_n_combinations(poly_dict, n-d, min_deg=d)
                for f in poly_dict[d]:
                    for g in g_list:
                        if f*g not in result:
                            result.append(f*g)
        return result

