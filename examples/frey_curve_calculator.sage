answer = dict()
for p in prime_range(6,20):
    DA = DiophantineAnalyzer(['a','b'])
    DA.add_restriction(CoprimeRestriction([a,b]))
    DA.add_restriction(CongruenceRestriction(b,[-1],4))
    f = a^p + b^p
    h = f.factor()[1][0]
    DA.add_restriction(PolynomialCongruenceRestriction(h,p^2,eq=False))
    DA.add_restriction(PolynomialPowerRestriction(h,5))
    L.<w> = CyclotomicField(p)
    answer[p] = dict()
    for Kdata in L.subfields():
        K = Kdata[0]
        S = K[a,b]
        fK = S(f)
        g_dict = make_degree_dict([d for (d,e) in fK.factor()])
        for g in find_degree_n_combinations(g_dict, 3):
            try:
                E = getEllipticCurveFromBinaryCubicForm(g)
                if K in answer[p]:
                    answer[p][K].append(E)
                else:
                    answer[p][K] = [E]
            except:
                pass
        for g in find_degree_n_combinations(g_dict, 6):
            try:
                E = getEllipticCurveFromBinaryCubicForm(g)
                if K in answer[p]:
                    answer[p][K].append(E)
                else:
                    answer[p][K] = [E]
            except:
                pass
