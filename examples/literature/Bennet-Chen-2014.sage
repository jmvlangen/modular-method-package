# Multi-Frey Q-curves and the Diophantine equation a^2 + b^6 = c^n
# Algebra & Number Theory volume 6 (2012), no.4

R.<a,b> = QQ[]
cp = a^2 + b^6

C = (CoprimeCondition([a,b]) &
     PowerCondition(cp, 5))

K.<i> = QuadraticField(-1)
L.<sqrtm3> = QuadraticField(-3)
G = K.galois_group()
E = FreyQcurve([0, 0, 0, -3*(5*b^3 + 4*a*i)*b, 2*(11*b^6 + 14*i*b^3*a - 2*a^2)],
               isogenies={G[0]: (QQ(1), 1), G[1]: (sqrtm3 , 3)},
               condition=C)

E = E.decomposable_twist()
print ""; print E.conductor()
# (4)*(1/4*izeta0^2 + 1)^n0*Rad_P( ((18195840*izeta0^3 - 145566720*izeta0 + 252129024)) * (b^3 + (-1/8*izeta0^3)*a) * (b^3 + (1/8*izeta0^3)*a)^3 )
#  where 
# n0 = 0 if ('a', 'b') is 1 of 24 possibilities mod 9
#      4 if ('a', 'b') is 1 of 48 possibilities mod 9
print ""; print E.conductor_restriction_of_scalars()
# 65536*3^(2*n0+4)*Norm(Rad_P( ((-72783360/7*izeta00zeta0^3 + 363916800/7*izeta00zeta0 + 252129024)) * (b^3 + (-1/28*izeta00zeta0^3 - 9/28*izeta00zeta0)*a) * (b^3 + (1/28*izeta00zeta0^3 + 9/28*izeta00zeta0)*a)^3 ))
#  where 
# n0 = 0 if ('a', 'b') is 1 of 24 possibilities mod 9
#      4 if ('a', 'b') is 1 of 48 possibilities mod 9
print ""; print E.newform_levels()
# [(48,)]  if ('a', 'b') is 1 of 24 possibilities mod 9
# [(432,)] if ('a', 'b') is 1 of 48 possibilities mod 9
