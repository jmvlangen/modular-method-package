#general setting
L = CyclotomicField(7)
R.<x,y> = L[]
F = x^7 + y^7

# Setting over totally real subfield
z7 = L.gens()[0]
K = L.subfield(z7 + z7^(-1),names = "alpha")[0]
Fa = K[x,y](F)
FaFactors = Fa.factor()

# The different degree 3 parts that we could look at
Fa1 = FaFactors[0][0] * FaFactors[1][0]
Fa2 = FaFactors[0][0] * FaFactors[2][0]
Fa3 = FaFactors[0][0] * FaFactors[3][0]
