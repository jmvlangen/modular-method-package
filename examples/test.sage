def randomCurve(K):
    R = K.maximal_order()
    try:
        return EllipticCurve([K(R.random_element()) for i in range(5)])
    except ArithmeticError:
        return randomCurve(K)
        
def same_local_data(d1, d2):
    return d1.bad_reduction_type() == d2.bad_reduction_type() \
       and d1.discriminant_valuation() == d2.discriminant_valuation() \
       and d1.tamagawa_number() == d2.tamagawa_number()

size = 1000
problems = []
for i in range(size):
    if (i+1) % 50 == 0:
        print i+1
    n = ZZ.random_element()
    while not n.is_squarefree() or n.is_square() or n == -6:
        n = ZZ.random_element()
    K = QuadraticField(n)
    #print K
    S = PolynomialRing(K, 1, names='v')
    ls = K.primes_above(2) + K.primes_above(3) + K.primes_above(5) + K.primes_above(7)
    p = ls[ZZ.random_element(len(ls))]
    #print p
    E = randomCurve(K)
    #print E
    cases = performTatesAlgorithm(E, prime=p, coefficient_ring=S, base_ring=K)
    answer = cases[0]
    if not same_local_data(answer, E.local_data(p)):
        problems.append((E, p))
    #print ""
print len(problems)
