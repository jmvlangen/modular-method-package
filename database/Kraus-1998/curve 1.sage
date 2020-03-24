# Article: Sur l'\'equation a^3 + b^3 = c^p
# Authors: Alain Kraus
# Journal: Experimental Mathematics, Volume 7, Issue 1, pages 1-13
# Source: https://www.tandfonline.com/doi/abs/10.1080/10586458.1998.10504355
#
# Equation: a^3 + b^3 = c^l
# Assumptions:
#   - a, b and c coprime
#   - l >= 17

# Conditions
R.<a, b> = QQ[]
cl = a^3 + b^3
C = (CoprimeCondition([a, b]) &
     PowerCondition(cl, 5)) # al = a^l

# The curve
# Y^2 = X^3 + 3*a*b*X + b^3 - a^3
a_invariants = [0, 0, 0, 3*a*b, b^3 - a^3]
E = FreyCurve(a_invariants, condition=C)

# Data
print("")
N = E.conductor()
print(N)
#
# (article: 2   * 3^2 * Rad if b = -1 (mod 4), c = 0 (mod 2)
#           2^4 * 3^2 * Rad if b =  1 (mod 4), c = 0 (mod 2)
#           2^4 * 3^2 * Rad if b = -1 (mod 4), c = 1 (mod 2)
#           2^3 * 3^2 * Rad if b =  1 (mod 4), c = 1 (mod 2), 2 || a
#           2^2 * 3^2 * Rad if b =  1 (mod 4), c = 1 (mod 2), 4 | a

# Newforms
print("")
nfs = E.newforms(algorithm='magma')
print(nfs)
# 
# (article: eliminated for specific l)
