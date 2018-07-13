Kcomp = Q.complete_field()
G = Kcomp.galois_group()
for beta in Q3.splitting_map('all'):
    flag = False
    for sigma in G:
        for tau in G:
            print (sigma, tau, beta(sigma) * beta(tau) * beta(sigma*tau)^(-1), Q3.c(sigma, tau))

### An overview of c
d = Q3.degree_map
for sigma in G:
    print [Q3.c(sigma, tau)/sqrt(d(sigma)*d(tau)*d(sigma*tau)^(-1)) for tau in G]

### Theta epsilon:
eps = Q3.splitting_character(galois=True)
for sigma in G:
    print [sqrt(eps(sigma)) * sqrt(eps(tau)) * sqrt(eps(sigma*tau))^(-1) for tau in G]

### Once more for good luck
Kcomp = Q1.complete_field()
Kl = Q1._Kl
G = Kcomp.galois_group()
d = Q1.degree_map
c = Q1.c
l = Q1.isogeny_lambda
eps = Q1.splitting_character(galois=True)
beta = Q1.splitting_map()
for sigma in G:
    print 'deps:', [sqrt(eps(sigma)) * sqrt(eps(tau)) * sqrt(eps(sigma*tau))^(-1) * sqrt(d(sigma) * d(tau) * d(sigma*tau)^(-1)) for tau in G]
    print '   c:', [c(sigma, tau) for tau in G]
    print '   l:', [l(sigma) * galois_field_change(sigma, Kl)(l(tau)) * l(sigma*tau)^(-1) for tau in G]
    print 'beta:', [beta(sigma) * beta(tau) * beta(sigma*tau)^(-1) for tau in G]
    print ""

### Sanity checks
ZNZ = [1,5,7,-1]
N = 12
d = Q1.degree_map
l = Q1.isogeny_lambda
beta = Q1.splitting_map()
c = Q1.c
L = CyclotomicField(N)
ls = [sqrt(L(a)) for a in [3,-3, -1]]
for n in ZNZ:
    s = cyclotomic_galois_isomorphism(n, N)
    print n, beta(s)
    if beta(s) not in QQ:
        print beta(s).minpoly()
    # for m in ZNZ:
    #     t = cyclotomic_galois_isomorphism(m, N)
    #     print (n,m), c(s,t)

### Comparing L-functions
Kcomp = Q.complete_field()
L = Q.splitting_image_field()
G = L.galois_group()
R.<x> = L[]
for p in prime_range(100):
    L_E = get_Lchi_E_p(E, p, x)
    print "p =", p
    print L_E
    for i in range(len(nfs2)):
        f = nfs2[i]
        ap = magma.Coefficient(f, p).sage().minpoly().change_ring(L).roots()[0][0]
        epsp = magma.DirichletCharacter(f)(p).sage()
        bp = p * epsp
        L_f = 1
        for sigma in G:
            L_f *= x^2 - sigma(ap) * x + bp
        if L_f == L_E:
            print i, L_f, "*"
        else:
            print i, L_f
    print ""

### Sanity checks:
f = nfs2[10]
for p in prime_range(5,500):
    L_E = get_Lchi_E_p(E, p, x)
    ap = magma.Coefficient(f, p).sage().minpoly().change_ring(L).roots()[0][0]
    epsp = magma.DirichletCharacter(f)(p).sage()
    bp = p * epsp
    L_f = 1
    for sigma in G:
        L_f *= x^2 - sigma(ap) * x + bp
    if L_f != L_E:
        print p
    print (1 - ap + bp).absolute_norm(), product(E.reduction(P).count_points() for P in K.primes_above(p))

### Finding the best possible twist
for p in prime_range(5, 10):
    M = magma.HeckeOperator(modsym, p) - (chi(p) * magma.Coefficient(f, p))
    K = magma.Kernel(M)
    print K

### My stuff
for f in nfs2:
    b = 0
    n = 0
    while b in ZZ:
        n += 1
        b = magma.Coefficient(f, n).sage()
    print b.parent().discriminant().squarefree_part()

### Sanity Check:
d = 

### Stuff
def prime_printer(a,b):
    ls = (a*b).prime_factors()
    if 2 not in ls:
        ls.append(2)
    for p in copy(ls):
        if hilbert_symbol(a,b,p) == 1:
            ls.remove(p)
    return ls

### More stuff
L = CyclotomicField(8)
sqrtm1 = sqrt(L(-1))
sqrt2 = sqrt(L(2))
sqrtm2 = sqrtm1 * sqrt2
ls = [a for a in range(24) if gcd(a,24) == 1]

def l(a):
    if a % 3 == 1:
        return 1
    else:
        return sqrtm2

def beta(a):
    if a % 12 == 1:
        return 1
    elif a % 12 == 5:
        return sqrtm2
    elif a % 12 == 7:
        return sqrtm1
    else:
        return sqrt2
    
def sl(a,b):
    if b % 3 == 1:
        return 1
    elif a % 8 == 1 or a % 8 == 3:
        return sqrtm2
    else:
        return -sqrtm2

for a in ls:
    print [l(a) * sl(a,b) * l(a*b % 24)^(-1) for b in ls]
print ""
for a in ls:
    print [beta(a) * beta(b) / beta(a*b % 24) for b in ls]
print ""
for a in ls:
    print [l(a) * sl(a,b) * l(a*b % 24)^(-1) / (beta(a) * beta(b) / beta(a*b % 24)) for b in ls]
