### Merel 1
# equation: a^n + b^n = 2*c^n
# source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
R.<ap, cp> = QQ[]
C = (CoprimeCondition([ap,cp]) &
     PowerCondition(ap, 7) & # ap = a^p
     PowerCondition(2*cp - ap, 7) & # bp = b^p
     PowerCondition(cp, 7) & # cp = c^p
     CongruenceCondition(ap + 1, 4)) # wlog
E = FreyCurve([0, (-ap)+(-2*cp), 0, (-ap)*(-2*cp), 0], # Y^2 = X*(X - ap)*(X - 2*cp)
              condition=C)
print ""
print E.conductor()
# 2^n0*Rad_P( 2^6 * cp^2 * ap^2 * (ap - 2*cp)^2 )
#  where 
# n0 = 5 if ('ap', 'cp') == (3, 1), (3, 3) mod 4
#      1 if ('ap', 'cp') is 1 of 8 possibilities mod 32
# (agrees with article)
print ""
print E.newforms()
# [(q - 2*q^5 + O(q^6), 'all')] if ('ap', 'cp') == (3, 1), (3, 3) mod 4
# (corresponds to the trivial solution (1, 1, 1))
# (agrees with article)

### Merel 2
# equation: a^n + b^n = c^2
# source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
# case: a*b even
R.<ap, c> = QQ[]
C = (CoprimeCondition([ap,c]) &
     PowerCondition(ap, 7) & # ap = a^p
     PowerCondition(c^2 - ap, 7) & # bp = b^p
     CongruenceCondition(ap, 2) & # wlog
     CongruenceCondition(c-1, 4)) # wlog
E = FreyCurve([1, (c - 1)/4, 0, ap / 2^6, 0], # Y^2 + X*Y = X^3 + (c-1)/4 * X^2 + a^p / 2^6 X
              condition=C)
print ""
print E.conductor()
# Rad_P( (1/4096) * ap^2 * (c^2 - ap) )
# (agrees with article)
print ""
print E.newforms()
# None
# (agrees with article)
del E

### Merel 3
# equation: a^n + b^n = c^2
# source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
# case: a*b odd
R.<ap, c> = QQ[]
C = (CoprimeCondition([ap,c]) &
     PowerCondition(ap, 7) & # ap = a^p
     PowerCondition(c^2 - ap, 7) & # bp = b^p
     CongruenceCondition(ap + 1, 4)) # wlog
E = FreyCurve([0, 2*c, 0, ap, 0], # Y^2 = X^3 + 2 * c * X^2 + a^p X
              condition=C)
print ""
print E.conductor()
# 32*Rad_P( 2^6 * ap^2 * (c^2 - ap) )
# (agrees with article)
print ""
print E.newforms()
# [(q - 2*q^5 + O(q^6), 'all')]
# (corresponds to the trivial solution (1, -1, 0))
# (agrees with article)
del E

### Merel 4
# equation: a^n + b^n = c^3
# source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
# case: c even
R.<bp, c> = QQ[]
C = (CoprimeCondition([bp,c]) &
     PowerCondition(c^3 - bp, 7) & # ap = a^p
     PowerCondition(bp, 7) & # bp = b^p
     CongruenceCondition(c, 2)) # wlog
E = FreyCurve([0, 0, bp, -3*((c/2)^3 + bp)*(c/2), -(c/2)^3*(2*(c/2)^3 - 5*bp)],
              # Y^2 + b^p*Y = X^3 - 3*((c/2)^3 + b^p)*(c/2)*X - (c/2)^3*(2*(c/2)^3 - 5*b^p)
              condition=C)
print ""
print E.conductor(verbose=2)
# 3^n0*Rad_P( (27) * bp * (c^3 - bp)^3 )
#  where 
# n0 = 2 if ('bp', 'c') is 1 of 710046 possibilities mod 2187
#      3 if ('bp', 'c') is 1 of 24 possibilities mod 9
#      1 if ('bp', 'c') is 1 of 1458 possibilities mod 2187
print ""
print E.newforms()
# [(q - 2*q^4 + O(q^6), 'all')] if ('bp', 'c') is 1 of 24 possibilities mod 9
# (corresponds to the trivial solution (1,-1,0))
# (agrees with article)
del E

### Merel 5
# equation: a^n + b^n = c^3
# source: http://www.math.mcgill.ca/darmon/pub/Articles/Research/18.Merel/pub18.pdf
# case: a*b even
R.<bp, c> = QQ[]
C = (CoprimeCondition([bp,c]) &
     PowerCondition(c^3 - bp, 7) & # ap = a^p
     PowerCondition(bp, 7) & # bp = b^p
     CongruenceCondition(c^3 - bp -1, 2) & # wlog a odd
     CongruenceCondition(bp, 2)) # wlog
E = FreyCurve([c, -c^2, 0, -3/2*c*bp, bp*((c^3 - bp) + 5/4*bp)],
              # Y^2 + c*X*Y = X^3 - c^2*X^2 - 3/2*c*b^p*X + b^p*(a^p + 5/4*b^p))
              condition=C)
print ""
print E.conductor(verbose=True)
# 3^n0*Rad_P( (27) * bp * (c^3 - bp)^3 )
#  where 
# n0 = 2 if ('bp', 'c') is 1 of 710046 possibilities mod 2187
#      3 if ('bp', 'c') is 1 of 24 possibilities mod 9
#      1 if ('bp', 'c') is 1 of 1458 possibilities mod 2187
print E.newforms()
# ???
del E

### Ellenberg 1
# equation: a^4 + b^2 = c^n
# source: http://www.math.wisc.edu/~ellenber/A4B2Cp.pdf
K.<i> = QuadraticField(-1)
R = K.ring_of_integers()
S.<a,b> = QQ[]
C = (CoprimeCondition([a,b]) &
     PowerCondition(a^4 + b^2, 3) & # a^4 + b^2 = c^p
     ~CongruenceCondition(b-1, 4)) # b != 1 mod 4
E = FreyCurve([0, 2*(1+i) * a, 0 , b + i*a^2, 0],
              parameter_ring=ZZ, condition=C)
# Y^2 = X^3 + 2*(1 + i)*a*X^2 + (b + i*a^2)*X
print ""
print E.conductor(verbose=True)
# Fractional ideal (i + 1)^n0*Rad_P( ((-64*i)) * (a^2 + (i)*b) * (a^2 + (-i)*b)^2 )
#  where 
# n0 = 12 if ('a', 'b') == (1, 0) mod 2
#      6  if ('a', 'b') == (0, 3), (2, 3) mod 4
# (agrees with article)
print ""
print E.newforms()
# ???
del E

