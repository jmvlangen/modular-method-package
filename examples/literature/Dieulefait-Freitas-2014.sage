# The fermat-type equations x^5 + y^5 = 2 z^p or 3 z^p solved through Q-curves.
# Mathematics of Computation, volume 83 (2014), number 286, pages 917-933

R.<a,b> = QQ[]
dcp = a^5 + b^5
phi = a + b
i = 0
while phi == a + b:
    phi = dcp.factor()[i][0]
    i += 1

K.<sqrt5> = QuadraticField(5)
(phi2, e2), (phi1, e1o) = phi.change_ring(K).factor()
w = (-1 + sqrt(5))/2
G.<sigma> = K.galois_group()

C = (CoprimeCondition([a,b]))

L.<sqrtm2> = QuadraticField(-2)
E = FreyQcurve([0, 2*(a + b), 0, - sigma(w)*phi1, 0],
               isogenies={G(1): (QQ(1), 1), sigma: (sqrtm2, 2)},
               condition=C)
E = E.decomposable_twist()

print ""; print E.decomposition_field()
# Number Field in zeta0 with defining polynomial x^4 - 5*x^2 + 5
print ""; print E.conductor()
# (2, -7/328*sqrt5sqrtm2zeta00^3 + 3/41*sqrt5sqrtm2zeta00^2 + 51/41*sqrt5sqrtm2zeta00 - 17/41)^n0*(5, -3/328*sqrt5sqrtm2zeta00^3 + 11/164*sqrt5sqrtm2zeta00^2 + 16/41*sqrt5sqrtm2zeta00 - 60/41)^n1*Rad_P( ((26800/41*sqrt5sqrtm2zeta00^3 + 132560/41*sqrt5sqrtm2zeta00^2 - 2302400/41*sqrt5sqrtm2zeta00 - 5224000/41)) * (a^2 + (1/82*sqrt5sqrtm2zeta00^3 - 1/164*sqrt5sqrtm2zeta00^2 - 35/41*sqrt5sqrtm2zeta00 - 43/41)*a*b + b^2) * (a^2 + (-1/82*sqrt5sqrtm2zeta00^3 + 1/164*sqrt5sqrtm2zeta00^2 + 35/41*sqrt5sqrtm2zeta00 + 2/41)*a*b + b^2)^2 )
#  where 
# n0 =  8 if ('a', 'b') is 1 of 6 possibilities mod 4  # a = 0 or b = 0 or a + b = 0 (mod 4)
#       6 if ('a', 'b') is 1 of 4 possibilities mod 4  # a = 2 or b = 2 (mod 4)
#       4 if ('a', 'b') is 1 of 4 possibilities mod 8  # a - b = 4 (mod 8)
#       0 if ('a', 'b') is 1 of 4 possibilities mod 8  # a - b = 0 (mod 8)
# n1 =  2 if ('a', 'b') is 1 of 20 possibilities mod 5 # a + b != 0 (mod 5)
#       0 if ('a', 'b') is 1 of 4 possibilities mod 5  # a + b = 0 (mod 5)
# NOTE: Does not agree with the literature
print ""; print E.conductor_restriction_of_scalars()
# 2^(2*n0+8)*5^(1*n1+6)*Norm(Rad_P( ((56960*zeta0^3 - 80960*zeta0^2 - 163200*zeta0 + 211200)) * (a^2 + (-zeta0^2 + 2)*a*b + b^2) * (a^2 + (zeta0^2 - 3)*a*b + b^2)^2 ))
#  where 
# n0 =  8 if ('a', 'b') is 1 of 6 possibilities mod 4  # a = 0 or b = 0 or a + b = 0 (mod 4)
#       6 if ('a', 'b') is 1 of 4 possibilities mod 4  # a = 2 or b = 2 (mod 4)
#       0 if ('a', 'b') is 1 of 4 possibilities mod 8  # a - b = 4 (mod 8)
#       4 if ('a', 'b') is 1 of 4 possibilities mod 8  # a - b = 0 (mod 8)
# n1 =  2 if ('a', 'b') is 1 of 20 possibilities mod 5 # a + b != 0 (mod 5)
#       0 if ('a', 'b') is 1 of 4 possibilities mod 5  # a + b = 0 (mod 5)
print ""; print E._newform_levels()
# [(1600, 1600)]             if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
# [(320, 1600), (1600, 320)] if ('a', 'b') is 1 of 6 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
# [(800, 800)]               if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 20 possibilities mod 5
# [(160, 800), (800, 160)]   if ('a', 'b') is 1 of 4 possibilities mod 4 and ('a', 'b') is 1 of 4 possibilities mod 5
# [(100, 100)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
# [(20, 100), (100, 20)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5
# [(400, 400)]               if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 20 possibilities mod 5
# [(80, 400), (400, 80)]     if ('a', 'b') is 1 of 4 possibilities mod 8 and ('a', 'b') is 1 of 4 possibilities mod 5
