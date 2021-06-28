### Preliminaries
from modular_method import *
from modular_method.elliptic_curves.twist import twist_elliptic_curve
from modular_method.diophantine_equations.conditions import ConditionalValue

### Function to compute in the Chen way
def chenMethod(E, nvars, P, nP, con):
    N = P^nP
    N = ZZ(N) if N in ZZ else ZZ(N.gens_reduced()[0])
    result = {}
    for vals in xmrange([N for i in range(nvars)]):
        if con(*vals):
            Evals = None
            while Evals is None:
                try:
                    Evals = E.specialize(vals)
                except ArithmeticError:
                    vals[0] += N
            e = Evals.local_data(P).conductor_valuation()
            if e in result:
                result[e].append(vals)
            else:
                result[e] = [vals]
    return result

### Function to compare the results
def same_answer(chenAnswer, frameworkAnswer):
    if isinstance(frameworkAnswer, ConditionalValue):
        if len(frameworkAnswer) != len(chenAnswer):
            return False
        for ans, con in frameworkAnswer:
            if ans not in chenAnswer:
                return False
            vals, N = con.pAdic_tree().give_as_congruence_condition()
            N = N.gens_reduced()[0]
            for val in chenAnswer[ans]:
                val = tuple(mod(x, N) for x in val)
                if val not in vals:
                    return False
        return True
    else:
        return (frameworkAnswer in chenAnswer) and (len(chenAnswer) == 1)

### Examples from Chen2010
R.<u, v> = QQ[]
s = v^2; t = u^2
K.<sqrt5> = QuadraticField(5)
delta = (-5 + 3*sqrt5)/2
kappa = (-1 + sqrt5)/2
nu = kappa^(-3)
coprime = CoprimeCondition([u, v])
Es = FreyQcurve([0, 0, 0,
                 -3*delta*((3 + 2*sqrt5)*s - 3*t),
                 4*v*((17 - 4*sqrt5)*s - (45 - 18*sqrt5)*t)],
                condition=coprime, guessed_degrees=[2])
Et = FreyQcurve([0, 0, 0,
                 -3*2^2*sqrt5*(3*s - (15 - 10*sqrt5)*t),
                 2^5*5*u*(9*s - (45 - 14*sqrt5)*t)],
                condition=coprime, guessed_degrees=[2])
S.<x> = K[]
Kbeta.<z> = K.extension(x^2 - (5 + sqrt5)/2)
gamma = 3*z^3 + 5*z^2 - 5*z - 5
Ebetas = Es.twist(gamma)
Ebetat = Es.twist(gamma)
print("Conductor of Ebetas from Chen2010 at P2")
P2 = Ebetas.definition_field().prime_above(2)
%time ans1 = chenMethod(Ebetas, 2, P2, 6, lambda u, v: not(2.divides(u) and 2.divides(v))) # 756 ms
%time ans2 = Ebetas.conductor_exponent(P2) # 521 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of Ebetas from Chen2010 at P3")
P3 = Ebetas.definition_field().prime_above(3)
%time ans1 = chenMethod(Ebetas, 2, P3, 4, lambda u, v: not (3.divides(u) and 3.divides(v))) # 57.7 s
%time ans2 = Ebetas.conductor_exponent(P3) # 792 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of Ebetas from Chen2010 at P5")
P5 = Ebetas.definition_field().prime_above(5)
zbetas = Ebetas.a4().factor().unit().factor()[1][0]
%time ans1 = chenMethod(FreyCurve(twist_elliptic_curve(Ebetas, zbetas^(-2))), 2, P5, 4, lambda u, v: not (5.divides(u) and 5.divides(v))) # 227 ms
%time ans2 = Ebetas.conductor_exponent(P5) # 110 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of Ebetat from Chen2010 at P2")
P2 = Ebetat.definition_field().prime_above(2)
%time ans1 = chenMethod(Ebetas, 2, P2, 6, lambda u, v: not(2.divides(u) and 2.divides(v))) # 668 ms
%time ans2 = Ebetat.conductor_exponent(P2) # 457 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of Ebetat from Chen2010 at P3")
P3 = Ebetat.definition_field().prime_above(3)
%time ans1 = chenMethod(Ebetat, 2, P3, 4, lambda u, v: not (3.divides(u) and 3.divides(v))) # 58.3 s
%time ans2 = Ebetat.conductor_exponent(P3) # 781 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of Ebetat from Chen2010 at P5")
P5 = Ebetat.definition_field().prime_above(5)
zbetat = Ebetat.a4().factor().unit().factor()[1][0]
%time ans1 = chenMethod(FreyCurve(twist_elliptic_curve(Ebetat, zbetat^(-2))), 2, P5, 4, lambda u, v: not (5.divides(u) and 5.divides(v))) # 216 ms
%time ans2 = Ebetat.conductor_exponent(P5) # 107 ms
print(same_answer(ans1, ans2)) # True

### Examples from Chen2012
R.<a, b> = QQ[]
K.<sqrt2> = QuadraticField(2)
coprime = CoprimeCondition([a, b])
E = FreyQcurve([0, 0, 0,
                -3*sqrt2*(4*a + 5*sqrt2*b^3)*b,
                4*sqrt2*(a^2 + 7*sqrt2*a*b^3 + 11*b^6)],
               condition=coprime, guessed_degrees=[3])
_.<sqrt6> = QuadraticField(6)
gamma = -3 + sqrt6
Ebeta = E.twist(gamma)
print("Conductor of Ebeta from Chen2012 at P2")
P2 = Ebeta.definition_field().prime_above(2)
%time ans1 = chenMethod(Ebeta, 2, P2, 8, lambda a, b: not (2.divides(a) and 2.divides(b))) # 497 ms
%time ans2 = Ebeta.conductor_exponent(P2) # 1.42 s
print(same_answer(ans1, ans2)) # True
print("Conductor of Ebeta from Chen2012 at P3")
P3 = Ebeta.definition_field().prime_above(3)
%time ans1 = chenMethod(Ebeta, 2, P3, 6, lambda a, b: not (3.divides(a) and 3.divides(b))) # 7.41 s
%time ans2 = Ebeta.conductor_exponent(P3) # 142 ms
print(same_answer(ans1, ans2)) # True

### Examples from BennettChen
R.<a, b> = QQ[]
K.<i> = QuadraticField(-1)
coprime = CoprimeCondition([a, b])
E = FreyQcurve([0, 0, 0,
                -3*(5*b^3 + 4*a*i)*b,
                2*(11*b^6 + 14*i*b^3*a - 2*a^2)],
               condition=coprime, guessed_degrees=[3])
_.<sqrtm3> = QuadraticField(-3)
gamma = (-3 + sqrtm3)/2
Ebeta = E.twist(gamma)
E_ = FreyCurve([0, 0, 0, 3*b^2, 2*a], condition=coprime)
print("Conductor of Ebeta from BennettChen at P2")
P2 = Ebeta.definition_field().prime_above(2)
%time ans1 = chenMethod(Ebeta, 2, P2, 8, lambda a, b: not (2.divides(a) and 2.divides(b))) # 4.14 s
%time ans2 = Ebeta.conductor_exponent(P2) # 6 min 25 s
print(same_answer(ans1, ans2)) # True
print("Conductor of Ebeta from BennettChen at P3")
P3 = Ebeta.definition_field().prime_above(3)
%time ans1 = chenMethod(Ebeta, 2, P3, 6, lambda a, b: not (3.divides(a) and 3.divides(b))) # 7.63 s
%time ans2 = Ebeta.conductor_exponent(P3) # 146 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E' from BennettChen at 2")
%time ans1 = chenMethod(E_, 2, 2, 8, lambda a, b: not (2.divides(a) and 2.divides(b))) # 18.1 s
%time ans2 = E_.conductor_exponent(2) # 13.9 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E' from BennettChen at 3")
%time ans1 = chenMethod(E_, 2, 3, 8, lambda a, b: not (3.divides(a) and 3.divides(b))) # 4h 1min 5s
%time ans2 = E_.conductor_exponent(3) # 137 ms
print(same_answer(ans1, ans2)) # True

### Examples from BennettChenDahmenYazdani
R.<s, t> = QQ[]
coprime = CoprimeCondition([s, t])
K.<sqrt3> = QuadraticField(3)
C1 = coprime & ~CongruenceCondition(s - t, 2) & ~CongruenceCondition(s - t, 3)
C2 = coprime & ~CongruenceCondition(s - t, 2) & ~CongruenceCondition(t, 3)
# For the case c odd
E1 = FreyCurve([0, 0, 0,
                 -6*(s^4 + 2*t*s^3 + 2*t^3*s + t^4),
                 -6*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)],
                condition=C1)
E2 = FreyQcurve([0, 2*(sqrt3 - 1)*(s - t), 0,
                  (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0],
                 condition=C1, guessed_degrees=[2])
print("Conductor of E1 (case c odd) from BennettChenYazdani at 2")
%time ans1 = chenMethod(E1, 2, 2, 6, lambda s, t: not (2.divides(s - t))) # 4.58 s
%time ans2 = E1.conductor_exponent(2) # 8.19 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E1 (case c odd) from BennettChenYazdani at 3")
%time ans1 = chenMethod(E1, 2, 3, 6, lambda s, t: not (3.divides(s - t))) # 2min 29s
%time ans2 = E1.conductor_exponent(3) # 9.74 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E2 (case c odd) from BennettChenYazdani at P2")
P2 = E2.definition_field().prime_above(2)
%time ans1 = chenMethod(E2, 2, P2, 6, lambda s, t: not (2.divides(s - t))) # 214 ms
%time ans2 = E2.conductor_exponent(P2) # 25.9 ms
print(same_answer(ans1, ans2)) # True
# For the case c even
E1a = FreyCurve([0, 0, 0,
                 -3*(-3*s^4 + 6*t^2*s^2 + t^4),
                 -12*s*t*(3*s^4 + t^4)],
                condition=C2)
E1b = FreyCurve([0, 0, 0,
                 -3*(3*s^4 + 6*t^2*s^2 - t^4),
                 -12*s*t*(3*s^4 + t^4)],
                condition=C2)
E2a = FreyQcurve([0, 4*(sqrt3 - 1)*t, 0,
                  -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 - sqrt3)*t^2), 0],
                 condition=C2, guessed_degrees=[2])
E2b = FreyQcurve([0, 4*(sqrt3 - 1)*t, 0,
                  -(sqrt3 - 1)^2*(sqrt3*s^2 + (-2 + sqrt3)*t^2), 0],
                 condition=C2, guessed_degrees=[2])
E3a = FreyQcurve([0, 12*(sqrt3 - 1)*s, 0,
                  3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 + 3)*s^2), 0],
                 condition=C2, guessed_degrees=[2])
E3b = FreyQcurve([0, 12*(sqrt3 - 1)*s, 0,
                  3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 - 3)*s^2), 0],
                 condition=C2, guessed_degrees=[2])
print("Conductor of E1a (case c even) from BennettChenYazdani at 2")
%time ans1 = chenMethod(E1a, 2, 2, 6, lambda s, t: not (2.divides(s - t))) # 818 ms
%time ans2 = E1a.conductor_exponent(2) # 19.9 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E1b (case c even) from BennettChenYazdani at 2")
%time ans1 = chenMethod(E1b, 2, 2, 6, lambda s, t: not (2.divides(s - t))) # 817 ms
%time ans2 = E1b.conductor_exponent(2) # 19.9 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E1a (case c even) from BennettChenYazdani at 3")
%time ans1 = chenMethod(E1a, 2, 3, 6, lambda s, t: not (3.divides(t))) # 2 min 27 s
%time ans2 = E1a.conductor_exponent(3) # 9.72 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E1b (case c even) from BennettChenYazdani at 3")
%time ans1 = chenMethod(E1b, 2, 3, 6, lambda s, t: not (3.divides(t))) # 2 min 29 s
%time ans2 = E1b.conductor_exponent(3) # 10.06 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E2a (case c even) from BennettChenYazdani at P2")
P2 = E2a.definition_field().prime_above(2)
%time ans1 = chenMethod(E2a, 2, P2, 6, lambda s, t: not (2.divides(s - t))) # 451 ms
%time ans2 = E2a.conductor_exponent(P2) # 40.7 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E2b (case c even) from BennettChenYazdani at P2")
P2 = E2b.definition_field().prime_above(2)
%time ans1 = chenMethod(E2b, 2, P2, 6, lambda s, t: not (2.divides(s - t))) # 383 ms
%time ans2 = E2b.conductor_exponent(P2) # 38.5 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E3a (case c even) from BennettChenYazdani at P2")
%time ans1 = chenMethod(E3a, 2, P2, 6, lambda s, t: not (2.divides(s - t))) # 382 ms
%time ans2 = E3a.conductor_exponent(P2) # 38.6 ms
print(same_answer(ans1, ans2)) # True
print("Conductor of E3b (case c even) from BennettChenYazdani at P2")
P2 = E3b.definition_field().prime_above(2)
%time ans1 = chenMethod(E3b, 2, P2, 6, lambda s, t: not (2.divides(s - t))) # 382 ms
%time ans2 = E3b.conductor_exponent(P2) # 39 ms
print(same_answer(ans1, ans2)) # True
