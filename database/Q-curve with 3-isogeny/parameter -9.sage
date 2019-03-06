load('~/Documents/SageFiles/load.sage')
load('~/Documents/SageFiles/examples/Q-curve families.sage')

print "Q-curve with 3-isogeny for parameter -9"

try:
    E = Qcurve_with_3_isogeny(-9).decomposable_twist()
    N = ZZ(E.conductor_restriction_of_scalars())
    primes = list(N.prime_factors()) + [-1]
    if 2 not in primes:
        primes.append(2)
    for p in primes:
        Etwist = E.twist(QQ(p))
        if Etwist.conductor_restriction_of_scalars() < N:
            E = Etwist
            N = E.conductor_restriction_of_scalars()

    print ""
    print "Decomposable curve of lowest conductor:"
    print ""
    print E
    print ""
    print "Newform levels:"
    print ""
    for tmp in E.newform_levels():
        print tmp
    print ""
    f, twists = E.newform()
    print "Newform of level %s:"%f.Level()
    print ""
    print f
    
except Exception as e:
    print ""
    print "An exception occured:"
    print e
