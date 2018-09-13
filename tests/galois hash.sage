keys = []
L = []
R = []
matches = []
for N in range(1, 100):
    for k in range(N):
        if gcd(k, N) == 1:
            Lv = hash((cyclotomic_galois_isomorphism(k, N=N), None))
            Rv = hash((cyclotomic_galois_isomorphism(k, N=N), N))
            if Lv in R:
                i = R.index(Lv)
                matches.append((k,N), keys[i])
            if Rv in L:
                i = L.index(Rv)
                matches.append(keys[i], (k, N))
            keys.append((k, N))
            L.append(Lv)
            R.append(Rv)
    print N, len(keys), len(matches)

###
@cached_function(key=lambda s, N: (0,(s, N)) if s in ZZ else (1,(str(s), s.parent().number_field(), N)))
def f(s, N=None):
    print "I got %s"%(s,)
    if not N is None:
        print "And I got %s"%(N,)

for N in range(1, 12):
    for k in range(N):
        if gcd(k, N) == 1:
            f(cyclotomic_galois_isomorphism(k, N=N), None)
            f(cyclotomic_galois_isomorphism(k, N=N), N)
