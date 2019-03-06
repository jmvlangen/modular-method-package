# Solving Fermat-type equations via modular Q-curves over polyquadratic fields
# ???

### Equation x^4 + 2*y^2 = z^p

R.<a,b> = QQ[]
cp = a^4 + 2*b^2

C = (CoprimeCondition([a,b]) &
     PowerCondition(cp, 5))

K.<r> = QuadraticField(-2)
G.<sigma> = K.galois_group()
E = FreyQcurve([0, 4*a, 0, 2*(a^2 + r*b), 0],
               isogenies={G(1): (QQ(1), 1), sigma: (r, 2)},
               condition=C)
print ""; print E.conductor()
# (r)^n0*Rad_P( (512) * (a^2 + (-r)*b) * (a^2 + (r)*b)^2 )
#  where 
# n0 = 12 if ('a', 'b') == (1, 1) mod 2
#      10 if ('a', 'b') == (1, 0) mod 2
print ""; print E.conductor_restriction_of_scalars()
# 2^(1*n0+6)*Norm(Rad_P( (512) * (a^2 + (-r)*b) * (a^2 + (r)*b)^2 ))
#  where 
# n0 = 12 if ('a', 'b') == (1, 1) mod 2
#      10 if ('a', 'b') == (1, 0) mod 2
print ""; print E.newform_levels()
# [(512,)] if ('a', 'b') == (1, 1) mod 2
# [(256,)] if ('a', 'b') == (1, 0) mod 2

### Equation x^4 + 3*y^2 = z^p

R.<a,b> = QQ[]
cp = a^4 + 3*b^2

C = (CoprimeCondition([a,b]) &
     PowerCondition(cp, 5))

K.<r> = QuadraticField(-3)
G.<sigma> = K.galois_group()
L.<sqrtm2> = QuadraticField(-2)
E = FreyQcurve([0, 4*a, 0, 2*(a^2 + r*b), 0],
               isogenies={G(1): (QQ(1), 1), sigma: (sqrtm2, 2)},
               condition=C)
E = E.decomposable_twist()
print ""; print E.conductor()
# (16)*Rad_P( ((-233610467125568/15*rsqrtm2rzeta0^7 + 583050224176000*rsqrtm2rzeta0^6 - 32212260985250176/15*rsqrtm2rzeta0^5 + 10494904035168000*rsqrtm2rzeta0^4 + 114728482001153792/15*rsqrtm2rzeta0^3 - 34983013450560000*rsqrtm2rzeta0^2 - 583609468822355456/15*rsqrtm2rzeta0 + 222330073930210304)) * (a^2 + (1/4800*rsqrtm2rzeta0^7 + 7/800*rsqrtm2rzeta0^5 + 63/400*rsqrtm2rzeta0^3 + 169/600*rsqrtm2rzeta0)*b) * (a^2 + (-1/4800*rsqrtm2rzeta0^7 - 7/800*rsqrtm2rzeta0^5 - 63/400*rsqrtm2rzeta0^3 - 169/600*rsqrtm2rzeta0)*b)^2 )
print ""; print E.conductor_restriction_of_scalars()
# 121029087867608368152576*Norm(Rad_P( ((-233610467125568/15*rsqrtm2rzeta0^7 + 583050224176000*rsqrtm2rzeta0^6 - 32212260985250176/15*rsqrtm2rzeta0^5 + 10494904035168000*rsqrtm2rzeta0^4 + 114728482001153792/15*rsqrtm2rzeta0^3 - 34983013450560000*rsqrtm2rzeta0^2 - 583609468822355456/15*rsqrtm2rzeta0 + 222330073930210304)) * (a^2 + (1/4800*rsqrtm2rzeta0^7 + 7/800*rsqrtm2rzeta0^5 + 63/400*rsqrtm2rzeta0^3 + 169/600*rsqrtm2rzeta0)*b) * (a^2 + (-1/4800*rsqrtm2rzeta0^7 - 7/800*rsqrtm2rzeta0^5 - 63/400*rsqrtm2rzeta0^3 - 169/600*rsqrtm2rzeta0)*b)^2 ))
print ""; print E.newform_levels()
# [(768, 768)]

# The levels of the article become visible over the cyclotomic field of order 48:

K = E.definition_field()
L = CyclotomicField(48)
iota = K.hom([a.minpoly().change_ring(L).roots()[0][0] for a in K.gens()], L)
Ebig = E.change_ring(iota)
print ""; print Ebig.conductor()
# (2, zeta48^2 + zeta48 + 1)^n0*Rad_P( ((-2462804146919424000*zeta48^14 - 1741465513021132800*zeta48^12 - 901448882318172160*zeta48^10 + 3364253029237596160*zeta48^6 + 3482931026042265600*zeta48^4 + 3364253029237596160*zeta48^2 + 3016306748181602304)) * (a^2 + (-2*zeta48^8 + 1)*b) * (a^2 + (2*zeta48^8 - 1)*b)^2 )
#  where 
# n0 = 12 if ('a', 'b') == (0, 1) mod 2
#      8  if ('a', 'b') == (1, 2), (3, 2) mod 4
#      0  if ('a', 'b') == (1, 0), (3, 0) mod 4
print ""; print Ebig.conductor_restriction_of_scalars()
# 2^(2*n0+96)*43046721*Norm(Rad_P( ((-2462804146919424000*zeta48^14 - 1741465513021132800*zeta48^12 - 901448882318172160*zeta48^10 + 3364253029237596160*zeta48^6 + 3482931026042265600*zeta48^4 + 3364253029237596160*zeta48^2 + 3016306748181602304)) * (a^2 + (-2*zeta48^8 + 1)*b) * (a^2 + (2*zeta48^8 - 1)*b)^2 ))
#  where 
# n0 = 12 if ('a', 'b') == (0, 1) mod 2
#      8  if ('a', 'b') == (1, 2), (3, 2) mod 4
#      0  if ('a', 'b') == (1, 0), (3, 0) mod 4
print ""; print Ebig.newform_levels()
# [(384, 768, 768, 384), (768, 384, 384, 768)]                                       if ('a', 'b') == (0, 1) mod 2
# [(192, 768, 768, 192), (768, 192, 192, 768)]                                       if ('a', 'b') == (1, 2), (3, 2) mod 4
# [(12, 768, 768, 192), (192, 768, 768, 12), (768, 24, 96, 768), (768, 96, 24, 768)] if ('a', 'b') == (1, 0), (3, 0) mod 4
