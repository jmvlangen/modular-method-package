S = PolynomialRing(QQ, ['x'])
S.inject_variables()
a = 9*x^4 - 92*x^3 - 42*x^2 - 60*x + 137
b = 101*x^6 + 30*x^5 + 795*x^4 - 2380*x^3 - 1605*x^2 + 654*x - 1627
E = EllipticCurve([0,0,0,-3*a,-2*b])
for p in prime_range(40):
    T = pAdicTree(list(S.gens()), p)
    T.apply_coprimality_restriction()
    T.apply_congruence_restriction(x, [6], 8)
    T.apply_congruence_anti_restriction(x, [1], 3)
    T.apply_congruence_anti_restriction(x, [19], 37)
    answers = performTatesAlgorithm(E, initial_values=T,
                                   only_calculate=['conductor','minimal_model'])
    print "Conductor valuations for %d:"%p
    for answer in answers:
        print answer[2].give_as_congruence_condition(), '->', answer[0]
        print 'Minimal model:', answer[1]
    print ""
    
# FLT:
# E = EllipticCurve([0,b-a,0,-a*b,0])
# T.apply_congruence_restriction(a, [-1], 4)
# T.apply_congruence_restriction(b, [0], 2^5)
#
# x^2 - 11 = y^l:
# E = EllipticCurve([0,-4*x,0,4*(x^2 - 11),0])
# T.apply_congruence_restriction(x, [0], 2)
# T.apply_congruence_anti_restriction(x, [0], 11)
#
# x^3 - x - 2 = y^l
# E = EllipticCurve([0,1,0,-x*(6 + x),-(2*x^3 + x^2 + 4*x + 4)])
# T.apply_congruence_restriction(x, [6], 8)
# 
# x^3 + 13 = y^l
# #(x even):
# E = EllipticCurve([0,0,0,-3*x^2,2*(x^3 + 26)])
# T.apply_congruence_restriction(x, [0], 2)
# #(x odd):
# E = EllipticCurve([0,0,0,-3*x^2,-2*(x^3 + 26)])
# T.apply_congruence_restriction(x, [11], 2^5)
# #(general):
# T.apply_congruence_anti_restriction(x, [-13], 3)
# T.apply_congruence_anti_restriction(x, [0], 13)
#
# x^4 + x^3 − 3x^2 + 11x + 2 = y^l:
# a = 9*x^4 - 92*x^3 - 42*x^2 - 60*x + 137
# b = 101*x^6 + 30*x^5 + 795*x^4 - 2380*x^3 - 1605*x^2 + 654*x - 1627
# E = EllipticCurve([0,0,0,-3*a,-2*b])
# T.apply_congruence_restriction(x, [6], 8)
# T.apply_congruence_anti_restriction(x, [1], 3)
# T.apply_congruence_anti_restriction(x, [19], 37)
#
