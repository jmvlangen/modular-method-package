### The (parametrized) curve
t = -3 # A non-square
K.<sqrt_t> = QuadraticField(t)
K2.<sqrt_m2> = QuadraticField(-2)
E = EllipticCurve([0,12,0,18*(sqrt_t + 1),0])
Q = Quer_invariants(E)
sigma = K.galois_group().gen()
Q.add_isogeny(sigma, 2, sqrt_m2)

### Doing the linear algebra
S = []
M, V, MT, VT, N, units = Q.c_beta_error('congruences', verbose=False, S=S)
M0 = minimal_echelon_form(M)
A = M0.right_kernel().basis_matrix().transpose()
M0T = MT * A
M1T = matrix([list(M0T[i]) + [VT[i]] for i in range(M0T.dimensions()[0])])
M2T = minimal_echelon_form(M1T, N=N)
if eliminate_zero_rows(M2T.delete_columns([M2T.dimensions()[1]-1])).dimensions()[0] < M2T.dimensions()[0]:
    print "Inconsistent system"
NI = N * matrix.identity(M2T.dimensions()[0])
M2 = M2T.change_ring(ZZ)
M3 = matrix([list(M2[i]) + list(NI[i]) for i in range(M2.dimensions()[0])])
C = M3.right_kernel().basis_matrix().transpose()
C0 = C[:M0T.dimensions()[1]]
C1 = C[M0T.dimensions()[1]:M2T.dimensions()[1]]
v0 = C1.solve_right(vector([-1]))
V0 = C1.right_kernel().basis_matrix().transpose()
v1 = A * C0 * v0
V1 = A * C0 * V0

### Constructing alpha, gamma and more
index, table = Q.c_beta_error('table', verbose=False)
alpha = {index[i] : product(units[j]^v1[i*len(units)+j] for j in range(len(units))) for i in range(len(index))}
x = 1
gamma = sum(s(x) / (alpha[s]^2) for s in index)

### A little searching script
index_dict = {v : k for k, v in enumerate(Q._c_beta_error_index)}
Kcomp = Q.complete_field()
gens = Kcomp.galois_group().gens()
orders = [g.order() for g in gens]

def alpha(s):
    if s not in a:
        i = max(i for i in range(len(s)) if s[i] != 0)
        s2 = [0]*len(s)
        s1 = list(s)
        s1[i] -= 1
        s2[i] += 1
        s1 = tuple(s1)
        s2 = tuple(s2)
        sigma1 = product(gens[j]^s1[j] for j in range(len(s)))
        i1 = index_dict[sigma1]
        i2 = index_dict[gens[i]]
        a[s] = alpha(s1) * sigma1(alpha(s2)) / QQ(Q._c_beta_error[i1][i2])
    return a[s]

def check_relation(s,t):
    sigma = product(gens[j]^s[j] for j in range(len(s)))
    i = index_dict[sigma]
    tau = product(gens[j]^t[j] for j in range(len(t)))
    j = index_dict[tau]
    st = tuple(s[j] + t[j] % orders[j] for j in range(len(s)))
    return alpha(s) * sigma(alpha(t)) == alpha(st) * QQ(Q._c_beta_error[i][j])

def my_iter(r, N):
    iters = [0]*r
    result = [0]*r
    for i in range(r):
        iters[i] = Kcomp.elements_of_bounded_height(N)
        while result[i] == 0 or abs(result[i].absolute_norm()) != 1:
            result[i] = iters[i].next()
    while True:
        yield result
        i = 0
        loop = True
        while loop:
            if i >= r:
                raise StopIteration
            try:
                result[i] = 0
                while result[i] == 0 or abs(result[i].absolute_norm()) != 1:
                    result[i] = iters[i].next()
                loop = False
            except StopIteration:
                iters[i] = Kcomp.elements_of_bounded_height(N)
                result[i] = 0
                while result[i] == 0 or abs(result[i].absolute_norm()) != 1:
                    result[i] = iters[i].next()
                i += 1

count = 0
for aprim in my_iter(len(gens), 5):
    a = dict()
    a[(0,0,0)] = QQ(Q._c_beta_error[0][0])
    a[(1,0,0)] = aprim[0]
    a[(0,1,0)] = aprim[1]
    a[(0,0,1)] = aprim[2]
    if product(a.itervalues()) != 0:
        count += 1
        for s in index_loop(orders):
            s = tuple(s)
            flag = True
            for t in [(1,0,0),(0,1,0),(0,0,1)]:
                flag = check_relation(s,t)
                if not flag:
                    break
            if not flag:
                break
        if flag:
            print "Solution found", aprim
        if count % 50 == 0:
            print "#", count

### New iterators:
def my_ZZ_iter(r, N):
    iters = [0]*r
    result = [0]*r
    for i in range(r):
        iters[i] = iter(ZZ.range(-N, N))
        result[i] = iters[i].next()
    while True:
        yield tuple(result)
        i = 0
        loop = True
        while loop:
            if i >= r:
                raise StopIteration
            try:
                result[i] = iters[i].next()
                loop = False
            except StopIteration:
                iters[i] = iter(ZZ.range(-N, N))
                result[i] = iters[i].next()
                i += 1

def norm_1_iterator(N):
    basis = Kcomp.integral_basis()
    for c in my_ZZ_iter(len(basis), N):
        x = sum(c[j] * basis[j] for j in range(len(basis)))
        if x != 0:
            try:
                n = x.absolute_norm().nth_root(8)
                x = x / n
                yield x
            except ValueError:
                pass

### Some counting:
count = 0
ls = []
for tmp in norm_1_iterator(10):
    count += 1
    ls.append(tmp)
    if count % 50 == 0:
        print "#", count
print "#", count
