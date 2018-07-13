DA = DiophantineAnalyzer(['ap','c']) # ap + bp = c
DA.add_restriction(CoprimeRestriction([ap,c]))
DA.add_restriction(PolynomialCongruenceRestriction(ap*(c^2-ap), 2))
DA.add_restriction(PolynomialPowerRestriction(ap, 7))
DA.add_restriction(PolynomialPowerRestriction(c^2-ap, 7))
DA.add_restriction(CongruenceRestriction(ap, 0, 2))
DA.add_restriction(CongruenceRestriction(c, 1, 4))
E = EllipticCurve([1, (c-1)/4, 0, ap/(2^6), 0])
answer = DA.compute_conductor(E)

###
DA = DiophantineAnalyzer(['bp','c'])
DA.add_restriction(CoprimeRestriction([bp,c]))
DA.add_restriction(PolynomialPowerRestriction(c^3 - bp, 7))
DA.add_restriction(PolynomialPowerRestriction(bp, 7))
DA.add_restriction(CongruenceRestriction(c, 0, 2))
E = EllipticCurve([0, 0, bp, -3/2 * ((c^3)/8 + bp) * c, -1/8 * c^3 * (1/4 * c^3 - 5 * bp)])
answer = DA.compute_conductor(E, model=True)
