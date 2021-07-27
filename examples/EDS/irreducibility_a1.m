// General setup
QQ := RationalField();
ZZ := Integers();
R<s, t> := PolynomialRing(QQ, 2);
S<a, b> := PolynomialRing(ZZ, 2);
A2 := AffinePlane(QQ);

// j-invariant of the Frey curve
jE := 2^6 * (3*t - 5)^3 / ((t - 1) * (t + 1)^2);
// where t = z / w^2 or -z / w^2

// Tries to find solutions to the Fermat equation
//     A*x^4 + B*y^4 + C*z^4 = 0
// by looking at the genus 1 quotient
//     A x^4 + B y^4 + C z^4 = 0 --> Y^2 Z = X^3 + (B C / A^2) X Z^2
//                   [x : y : z] |-> [B y^2 z : B x^2 y : -A z^3 ]
// If the quotient curve has rank 0 it will determine all triples [x, y, z]
// that map to its torsion points.
// The return values in that case are true and the list of those triples.
// Otherwise the return values are false and undefined.
function signature444quotientSolver(A, B, C)
    E := EllipticCurve([0, 0, 0, B*C / A^2, 0]);
    if Rank(E) gt 0 then
        return false, _;
    end if;
    Etor, map := TorsionSubgroup(E);
    // For a point P := (a : b : c)
    // Find the triples of coprime integers x, y, z that solve
    //    B y^2 z = a k
    //    B x^2 y = b k
    //     -A z^3 = c k
    // for some non-zero k. Note that these x, y, z satisfy
    //            (y/z)^2 = (- a A) / (B c)
    //    (x/z)^2 (y/z)   = (- b A) / (B c)
    // so we find x, y, z by first solving the latter two over Q
    result := [];
    for Ptor in Etor do
        P := map(Ptor);
        if P[3] eq 0 then
            Append(~result, [1, 1, 0]);
            Append(~result, [-1, 1, 0]);
            Append(~result, [1, -1, 0]);
            Append(~result, [-1, -1, 0]);
        else
            ydivzsq := - P[1] * A / (B * P[3]);
            test, ydivz := IsSquare(ydivzsq);
            if test then
                if ydivz eq 0 then
                    Append(~result, [1, 0, 1]);
                    Append(~result, [-1, 0, 1]);
                    Append(~result, [1, 0, -1]);
                    Append(~result, [-1, 0, -1]);
                else
                    xdivzsq := (- P[2]*A) / (B * P[3] * ydivzsq);
                    test, xdivz := IsSquare(xdivzsq);
                    if test then
                        if xdivz eq 0 then
                            Append(~result, [0, 1, 1]);
                            Append(~result, [0, -1, 1]);
                            Append(~result, [0, 1, -1]);
                            Append(~result, [0, -1, -1]);
                        else
                            z := Lcm(Denominator(xdivz), Denominator(ydivz));
                            y := Integers() ! (ydivz * z);
                            x := Integers() ! (xdivz * z);
                            Append(~result, [ x,  y,  z]);
                            Append(~result, [ x,  y, -z]);
                            Append(~result, [ x, -y,  z]);
                            Append(~result, [ x, -y, -z]);
                            Append(~result, [-x,  y,  z]);
                            Append(~result, [-x,  y, -z]);
                            Append(~result, [-x, -y,  z]);
                            Append(~result, [-x, -y, -z]);
                        end if;
                    end if;
                end if;
            end if;
        end if;
    end for;
    return true, {xyz : xyz in result};
end function;

// Tries to find solutions to the Fermat equation
//     A*x^4 + B*y^4 + C*z^4 = 0
// by the following steps:
// First determine if the associated curve is locally solvable at 2, 3 and 5
// If not returns true and an empty list
// Otherwise it will look at the different possible genus 1 quotients
//     A x^4 + B y^4 + C z^4 = 0 --> Y^2 Z = X^3 + (B C / A^2) X Z^2
//                   [x : y : z] |-> [B y^2 z : B x^2 y : -A z^3 ]
// If at least one of them has rank 0, will determine all the torsion
// points thereof.
// In that case this function returns True and a list of triples [x, y, z]
// that correspond to the points mapping to these torsion points.
// Otherwise this function returns False and the second return value is undefined.
function signature444solver(A, B, C)
    R<x, y, z> := PolynomialRing(QQ, 3);
    P2 := ProjectivePlane(QQ);
    curve := Curve(P2, A*x^4 + B*y^4 + C*z^4);
    for p in [2, 3, 5] do
        if not IsLocallySolvable(curve, p) then
            return true, [];
        end if;
    end for;
    result, list := signature444quotientSolver(A, B, C);
    if result then
        return result, list;
    end if;
    result, list := signature444quotientSolver(A, C, B);
    if result then
        return result, [[P[1], P[3], P[2]] : P in list];
    end if;
    result, list := signature444quotientSolver(B, A, C);
    if result then
        return result, [[P[2], P[1], P[3]] : P in list];
    end if;
    result, list := signature444quotientSolver(B, C, A);
    if result then
        return result, [[P[2], P[3], P[1]] : P in list];
    end if;
    result, list := signature444quotientSolver(C, A, B);
    if result then
        return result, [[P[3], P[1], P[2]] : P in list];
    end if;
    result, list := signature444quotientSolver(C, B, A);
    if result then
        return result, [[P[3], P[2], P[1]] : P in list];
    end if;
    return false, _;
end function;

// Check that in case we have a reducible representation at 3
// then we find that
//    c^2 D B^4 = 2^8 * a^3 * b^3 * (a - b) * (a + 8*b)
// for some coprime a, b in \Z
// and some c in \Z only divisible by prime numbers 2 and 3
// Note that this part is independent of the chosen D!
X03 := SmallModularCurve(3);
assert Genus(X03) eq 0;
j3 := jInvariant(X03, 3);
j3 := Evaluate(j3, [s, 1]);
f3 := Numerator(j3) * Denominator(jE) - Numerator(jE) * Denominator(j3);
C3 := Curve(A2, f3);
C3 := ProjectiveClosure(C3);
assert Genus(C3) eq 0;
P3 := C3 ! [0, 1];
phi3 := Parametrization(C3, P3);
t3 := DefiningPolynomials(phi3)[2] / DefiningPolynomials(phi3)[3];
z31 := 2187*a^2 - 230688*a*b + 6079232*b^2;
z3sq := 3 * z31^2;
w31 := 2187*a^2 - 234144*a*b + 6269696*b^2;
w32 := 6561*a^2 - 681696*a*b + 17682688*b^2;
w3 := w31 * w32;
assert t3 eq w3 / z3sq;
// We replace with a better model
abnew := [32*a + 4336*b, 81*b];
z3 := a^2 + 4*a*b - 8*b^2;
w31 := a^2 + 8*b^2;
w32 := a^2 + 8*a*b - 8*b^2;
w3 := w31 * w32;
t3 := Evaluate(t3, abnew);
assert t3 eq w3 / z3^2;
// Note that we must have some coprime integers a and b
// and a non-zero integer c such that
// c z^2 = z3(a, b)^2
// c w = w31(a, b) * w32(a, b)
// So c divides the resultant of z3^2 and w3
assert Evaluate(Resultant(z3^2, w3, a), [1, 1]) eq 2^28 * 3^2;
// Finally this shows the equation we wanted to show
assert w3^2 - z3^4 eq 2^8 * a^3 * b^3 * (a - b) * (a + 8*b);

// A function that given a specific D
// computes whether the equation
//    c^2 D B^4 = 2^8 * a^3 * b^3 * (a - b) * (a + 8*b) [*3]
// for some coprime a, b in \Z
// and some c in \Z only divisible by prime numbers 2 and 3
// does not have solutions
// which by the previous implies the mod 3 representation is irreducible
function irreducible_at_3(D)
    // Note that equation [*3] implies that
    // a, b, a - b, and a + 8*b are constants times fourth powers
    // Below we compute the possible constants a_1, ..., a_4 s.t.
    //      a       = a_1 c_1^4
    //            b = a_2 c_2^4
    //      a -   b = a_3 c_3^4
    //      a + 8*b = a_4 c_4^4
    // We first compute the possibilities for powers of (-1) and
    // each prime p appearing in these constants
    // Note that 2 and 3 are special, as these are the primes
    // that might divide c. Furthermore we have
    //    gcd(a, b) = gcd(a, a - b) = gcd(b, a - b) = gcd(b, a + 8*b) = 1
    //    gcd(a, a + 8*b) divides 8
    //    gcd(a - b, a + 8*b) divides 9
    // so 2 and 3 are the only primes that can appear in two of the factors
    // a, b, a - b and a + 8*b simultaneously.
    cfs := [{[(-1)^a1, (-1)^a2, (-1)^a3, (-1)^a4] :
	         a1 in [0..1], a2 in [0..1], a3 in [0..1], a4 in [0..1] |
	         (a1 + a2 + a3 + a4 - ((D lt 0) select 1 else 0)) mod 2 eq 0}];
    for p in PrimeDivisors(2*3*D) do
        if p eq 2 then
	        // ord_2(z3(a, b))  = 0 if ord_2(a) = 0
	        //                    2 if ord_2(a) = 1
	        //                    3 if ord_2(a) > 1
	        // ord_2(w31(a, b)) = 0 if ord_2(a) = 0
	        //                    2 if ord_2(a) = 1
	        //                    3 if ord_2(a) > 1
	        // ord_2(w32(a, b)) = 0 if ord_2(a) = 0
	        //                    2 if ord_2(a) = 1
	        //                    3 if ord_2(a) > 1
	        // Therefore
	        // ord_2(c)         = 0 if ord_2(a) = 0
	        // ord_2(c)         = 4 if ord_2(a) = 1
	        // ord_2(c)         = 6 if ord_2(a) > 1
	        // Furthermore we have the relation
	        // 2 ord_2(c) + ord_2(D) = 3 ord_2(a) + 3 ord_2(b) + ord_2(a - b) + ord_2(a + 8*b) mod 4
	        // hence
	        // ord_2(D) = 3 ord_2(a1) + 3 ord_2(a2) + ord_2(a3) + ord_2(a4) mod 4
	        Append(~cfs, {[p^c1, p^c2, p^c3, p^c4] :
		                  c1 in [0..3], c2 in [0..3], c3 in [0..3], c4 in [0..3] |
		                  ((c2 + c3 eq 0) or (c1 + c3 + c4 eq 0) or (c1 + c2 + c4 eq 0)) and
		                  ((3*c1 + 3*c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0)});
        elif p eq 3 then
	        // Note z3(a, b) = (a + 2*b)^2 - 12*b^2, so
	        // ord_3(z3(a, b))  = 0 if ord_3(a - b) = 0
	        //                    1 if ord_3(a - b) > 0
	        // Note w31(a, b) = (a - b)*(a + b) mod 9, so
	        // ord_3(w31(a, b)) = 0 if ord_3(a - b) = 0 and ord_3(a + b) = 0
	        //                  > 0 otherwise
	        // Note w32(a, b) = (a + 4*b)^2 - 24*b^2, so
	        // ord_3(w32(a, b)) = 0 if ord_3(a + b) = 0
	        //                    1 if ord_3(a + b) > 0
	        // Therefore
	        // ord_3(c)         = 0 if ord_2(a - b) = 0
	        //                  = 2 if ord_2(a - b) > 0
	        // Furthermore we have the relation
	        // 2 ord_3(c) + ord_3(D) = 3 ord_3(a) + 3 ord_3(b) + ord_3(a - b) + ord_3(a + 8*b) mod 4
	        // hence
	        // ord_3(D) = 3 ord_3(a1) + 3 ord_3(a2) + ord_3(a3) + ord_3(a4) mod 4
	        Append(~cfs, {[p^c1, p^c2, p^c3, p^c4] :
		                  c1 in [0..3], c2 in [0..3], c3 in [0..3], c4 in [0..3] |
		                  ((c2 + c3 + c4 eq 0) or (c1 + c3 + c4 eq 0) or (c1 + c2 eq 0)) and
		                  ((3*c1 + 3*c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0)});
        else
	        // We have ord_p(c) = 0, so just the relation
	        // ord_p(D) = 3 ord_p(a1) + 3 ord_p(a2) + ord_p(a3) + ord_p(a4) mod 4
	        Append(~cfs, {[p^c1, p^c2, p^c3, p^c4] :
		                  c1 in [0..3], c2 in [0..3], c3 in [0..3], c4 in [0..3] |
		                  ((c2 + c3 + c4 eq 0) or (c1 + c3 + c4 eq 0) or
		                   (c1 + c2 + c4 eq 0) or (c1 + c2 + c3 eq 0)) and
		                  ((3*c1 + 3*c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0)});
        end if;
    end for;
    cfs := CartesianProduct(cfs);
    cfs := [ [ &* [ cf[i] : cf in cfls] : i in [1..4] ] : cfls in cfs];
    // Using the linear relations between a, b, a - b, and a + 8*b
    // we can construct four curves for each possibility of the
    // form a*x^4 + b*y^4 = c*z^4
    // These relations are
    //       b      +   (a - b)   =    a
    //       a      + 8* b        =   (a + 8*b)
    //    8*(a - b) +   (a + 8*b) = 9* a
    //    9* b      +   (a - b)   =   (a + 8*b)
    bad_cfs := {};
    bad_ab := {};
    for cf in cfs do
        test, xyzls := signature444solver(cf[2], cf[3], -cf[1]);
        if test then
            for xyz in xyzls do
                a := cf[1]*xyz[3]^4;
                b := cf[2]*xyz[1]^4;
                if a - b eq cf[3]*xyz[2]^4 then
                    Include(~bad_ab, [a, b]);
                end if;
            end for;
        else
            test, xyzls := signature444solver(cf[1], 8*cf[2], -cf[4]);
            if test then
                for xyz in xyzls do
                    a := cf[1]*xyz[1]^4;
                    b := cf[2]*xyz[2]^4;
                    if a + 8*b eq cf[4]*xyz[3]^4 then
                        Include(~bad_ab, [a, b]);
                    end if;
                end for;
            else
                test, xyzls := signature444solver(8*cf[3], cf[4], -9*cf[1]);
                if test then
                    for xyz in xyzls do
                        a := cf[1]*xyz[3]^4;
                        b := a - cf[3]*xyz[1]^4;
                        if a + 8*b eq cf[4]*xyz[2]^4 then
                            Include(~bad_ab, [a, b]);
                        end if;
                    end for;
                else
                    test, xyzls := signature444solver(9*cf[2], cf[3], -cf[4]);
                    if test then
                        for xyz in xyzls do
                            b := cf[2]*xyz[1]^4;
                            a := b + cf[3]*xyz[2]^4;
                            if a + 8*b eq cf[4]*xyz[3]^4 then
                                Include(~bad_ab, [a, b]);
                            end if;
                        end for;
                    else
                        Include(~bad_cfs, cf);
                    end if;
                end if;
            end if;
        end if;
    end for;
    if # bad_cfs gt 0 then
        return false, [*bad_cfs, bad_ab*];
    end if;
    if # bad_ab eq 0 then
        return true, _;
    end if;
    bad_zwB := {};
    for ab in bad_ab do
        a := ab[1];
        b := ab[2];
        czsq := (a^2 + 4*a*b - 8*b^2)^2;
        cw := (a^2 + 8*b^2) * (a^2 + 8*a*b - 8*b^2);
        c := Gcd(czsq, cw);
        zsq := Integers() ! (czsq / c);
        w := Integers() ! (cw / c);
        test, z := IsSquare(zsq);
        if test then
            DB4 := w^2 - z^4;
            if DB4 mod D eq 0 then
                B4 := Integers() ! (DB4 / D);
                test, B2 := IsSquare(B4);
                if test then
                    test, B := IsSquare(B2);
                    if test then
                        Include(~bad_zwB, [z, w, B]);
                    end if;
                end if;
            end if;
        end if;
    end for;
    ED := EllipticCurve([0, 0, 0, D, 0]);
    bad_points := {ED ! [zwB[1]^2*zwB[3], zwB[1]*zwB[2], zwB[3]^3] : zwB in bad_zwB};
    if bad_points eq { ED ! [0, 1, 0] } then
        return true, bad_points;
    else
        return false, bad_points;
    end if;
end function;

// We apply the function reducible_at_3 to various values of D
irreducible_at_3(-2);
// true { (0 : 1 : 0) }
irreducible_at_3(3);
// true { (0 : 1 : 0) }
irreducible_at_3(-17);
// true { (0 : 1 : 0) }
irreducible_at_3(125);
// false { (121/4 : 1419/8 : 1), (0 : 1 : 0), (121/4 : -1419/8 : 1) }

// Check that in case we have a reducible representation at 5
// then we find that
//    c^2 D B^4 = 2^8 * a^5 * b^5 * (a - b) * (a + 4*b)
// for some coprime a, b in \Z
// and some c in \Z only divisible by prime numbers 2 and 5
// Note that this part is independent of the chosen D!
X05 := SmallModularCurve(5);
assert Genus(X05) eq 0;
j5 := jInvariant(X05, 5);
j5 := Evaluate(j5, [s, 1]);
f5 := Numerator(j5) * Denominator(jE) - Numerator(jE) * Denominator(j5);
C5 := Curve(A2, f5);
C5 := ProjectiveClosure(C5);
assert Genus(C5) eq 0;
P5 := C5 ! [0, 1];
phi5 := Parametrization(C5, P5);
t5 := DefiningPolynomials(phi5)[2] / DefiningPolynomials(phi5)[3];
z51 := 1625625*a^2 - 138748886400*a*b + 2960592618409984*b^2;
z52 := 1625625*a^2 - 138746438400*a*b + 2960488141287424*b^2;
z5sq := z51 * z52^2;
w5 := 4295968701416015625*a^6 - 1099984059625938750000000*a^5*b +
      117354685225489694136000000000*a^4*b^2 -
      6677489319085419352663080960000000*a^3*b^3 +
      213721639890166222651532224128614400000*a^2*b^4 -
      3648232615051429955188489068117958538035200*a*b^5 +
      25948084500870947622494893023137858568620867584*b^6;
assert t5 eq w5 / z5sq;
// We replace with a better model
abnew := [960*a + 54411328*b, 1275*b];
z51 := a^2 + 4*b^2;
z52 := a^2 + 2*a*b - 4*b^2;
z5sq := z51 * z52^2;
w5 := a^6 + 4*a^5*b + 64*a*b^5 - 64*b^6;
t5 := Evaluate(t5, abnew);
assert t5 eq w5 / z5sq;
// Note that we must have some coprime integers a and b
// and a non-zero integer c such that
// c z^2 = z51(a, b) * z52(a, b)^2
// c w   = w5(a, b)
// So c divides the resultant of z5sq and w5
assert Evaluate(Resultant(z5sq, w5, a), [1, 1]) eq 2^42 * 5;
// Finally this shows the equation we wanted to show
assert w5^2 - z5sq^2 eq 2^8 * a^5 * b^5 * (a - b) * (a + 4*b);

// A function that given a specific D
// computes whether the equation
//    c^2 D B^4 = 2^8 * a^5 * b^5 * (a - b) * (a + 4*b) [*5]
// for some coprime a, b in \Z
// and some c in \Z only divisible by prime numbers 2 and 5
// does not have solutions
// which by the previous implies the mod 5 representation is irreducible
function irreducible_at_5(D)
    // Note that equation [*5] implies that
    // a, b, a - b, and a + 4*b are constants times fourth powers
    // Below we compute the possible constants a_1, ..., a_4 s.t.
    //      a       = a_1 c_1^4
    //            b = a_2 c_2^4
    //      a -   b = a_3 c_3^4
    //      a + 4*b = a_4 c_4^4
    // We first compute the possibilities for powers of (-1) and
    // each prime p appearing in these constants
    // Note that 2 and 5 are special, as at these are the primes
    // that might divide c. Furthermore we have
    //    gcd(a, b) = gcd(a, a - b) = gcd(b, a - b) = gcd(b, a + 4*b) = 1
    //    gcd(a, a + 4*b) divides 4
    //    gcd(a - b, a + 4*b) divides 5
    // so 2 and 5 are the only primes that can appear in two of the factors
    // a, b, a - b and a + 4*b simultaneously.
    cfs := [{[(-1)^a1, (-1)^a2, (-1)^a3, (-1)^a4] :
	         a1 in [0..1], a2 in [0..1], a3 in [0..1], a4 in [0..1] |
	         (a1 + a2 + a3 + a4 - ((D lt 0) select 1 else 0)) mod 2 eq 0}];
    for p in PrimeDivisors(2*5*D) do
        if p eq 2 then
	        // ord_2(z51(a, b)) = 0 if ord_2(a) = 0
	        //                    3 if ord_2(a) = 1
	        //                    2 if ord_2(a) > 1
	        // Note z52(a, b) = (a + b)^2 - 5*b^2, so
	        // ord_2(z52(a, b)) = 0 if ord_2(a) = 0
	        //                    2 if ord_2(a) > 0
	        // ord_2(w5(a, b))  = 0 if ord_2(a) = 0
	        //                 >= 7 if ord_2(a) = 1
	        //                  = 6 if ord_2(a) > 1
	        // Therefore
	        //         ord_2(c) = 0 if ord_2(a) = 0
	        //                    7 if ord_2(a) = 1
	        //                    6 if ord_2(a) > 1
	        // Furthermore we have the relation
	        // 2 ord_2(c) + ord_2(D) = 5 ord_2(a) + 5 ord_2(b) + ord_2(a - b) + ord_2(a + 4*b) mod 4
	        // hence if ord_2(a) = 1
	        //   2 + ord_2(D) = ord_2(a1) + ord_2(a2) + ord_2(a3) + ord_2(a4) mod 4
	        // and otherwise
	        //   ord_2(D) = ord_2(a1) + ord_2(a2) + ord_2(a3) + ord_2(a4) mod 4
	        Append(~cfs, {[p^c1, p^c2, p^c3, p^c4] :
		                  c1 in [0..3], c2 in [0..3], c3 in [0..3], c4 in [0..3] |
		                  ((c2 + c3 eq 0) or (c1 + c3 + c4 eq 0) or (c1 + c2 + c4 eq 0)) and
		                  ((((c1 eq 1) and ((2 + c1 + c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0)) or
			                ((c1 + c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0)))});
        elif p eq 5 then
	        // Note z51(a, b) = (a - b)*(a + b) mod 5, so
	        // ord_5(z51(a, b)) = 0 if ord_5(a - b) = 0 and ord_5(a + b) = 0
	        //                  > 0 otherwise
	        // Note z52(a, b) = (a + b)^2 - 5*b^2, so
	        // ord_5(z52(a, b)) = 0 if ord_5(a + b) = 0
	        //                    1 if ord_5(a + b) > 0
	        // Note w5(a, b) = (a - b)^6 mod 5, so
	        // ord_5(w5(a, b))  = 0 if ord_5(a - b) = 0
	        //                  > 0 if ord_5(a - b) > 0
	        // Note that the resultant of z5sq and w5 only had one factor 5.
	        // Therefore
	        //         ord_5(c) = 0 if ord_5(a - b) = 0
	        //                    1 if ord_5(a - b) > 0
	        // Furthermore we have the relation
	        // 2 ord_5(c) + ord_5(D) = 5 ord_5(a) + 5 ord_5(b) + ord_5(a - b) + ord_5(a + 4*b) mod 4
	        // hence if ord_5(a - b) > 0
	        //   2 + ord_5(D) = ord_5(a1) + ord_5(a2) + ord_5(a3) + ord_5(a4) mod 4
	        // and otherwise
	        //   ord_5(D) = ord_5(a1) + ord_5(a2) + ord_5(a3) + ord_5(a4) mod 4
	        Append(~cfs, {[p^c1, p^c2, p^c3, p^c4] :
		                  c1 in [0..3], c2 in [0..3], c3 in [0..3], c4 in [0..3] |
		                  ((c2 + c3 + c4 eq 0) or (c1 + c3 + c4 eq 0) or (c1 + c2 eq 0)) and
		                  ((((c1 + c2 eq 0) and ((2 + c1 + c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0)) or
			                ((c3 eq 0) and ((c1 + c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0))))});
        else
	        // We have ord_p(c) = 0, so just the relation
	        // ord_p(D) = ord_p(a1) + ord_p(a2) + ord_p(a3) + ord_p(a4) mod 4
	        Append(~cfs, {[p^c1, p^c2, p^c3, p^c4] :
		                  c1 in [0..3], c2 in [0..3], c3 in [0..3], c4 in [0..3] |
		                  ((c2 + c3 + c4 eq 0) or (c1 + c3 + c4 eq 0) or
		                   (c1 + c2 + c4 eq 0) or (c1 + c2 + c3 eq 0)) and
		                  ((c1 + c2 + c3 + c4 - Valuation(D, p)) mod 4 eq 0)});
        end if;
    end for;
    cfs := CartesianProduct(cfs);
    cfs := [ [ &* [ cf[i] : cf in cfls] : i in [1..4] ] : cfls in cfs];
    // Using the linear relations between a, b, a - b, and a + 4*b
    // we can construct four curves for each possibility of the
    // form a*x^4 + b*y^4 = c*z^4
    // These relations are
    //      b      +   (a - b)   =    a
    //      a      + 4* b        =   (a + 4*b)
    //   4*(a - b) +   (a + 4*b) = 5* a
    //   5* b      +   (a - b)   =   (a + 4*b)
    bad_cfs := {};
    bad_ab := {};
    for cf in cfs do
        test, xyzls := signature444solver(cf[2], cf[3], -cf[1]);
        if test then
            for xyz in xyzls do
                a := cf[1]*xyz[3]^4;
                b := cf[2]*xyz[1]^4;
                if a - b eq cf[3]*xyz[2]^4 then
                    Include(~bad_ab, [a, b]);
                end if;
            end for;
        else
            test, xyzls := signature444solver(cf[1], 4*cf[2], -cf[4]);
            if test then
                for xyz in xyzls do
                    a := cf[1]*xyz[1]^4;
                    b := cf[2]*xyz[2]^4;
                    if a + 4*b eq cf[4]*xyz[3]^4 then
                        Include(~bad_ab, [a, b]);
                    end if;
                end for;
            else
                test, xyzls := signature444solver(4*cf[3], cf[4], -5*cf[1]);
                if test then
                    for xyz in xyzls do
                        a := cf[1]*xyz[3]^4;
                        b := a - cf[3]*xyz[1]^4;
                        if a + 4*b eq cf[4]*xyz[2]^4 then
                            Include(~bad_ab, [a, b]);
                        end if;
                    end for;
                else
                    test, xyzls := signature444solver(5*cf[2], cf[3], -cf[4]);
                    if test then
                        for xyz in xyzls do
                            b := cf[2]*xyz[1]^4;
                            a := b + cf[3]*xyz[2]^4;
                            if a + 4*b eq cf[4]*xyz[3]^4 then
                                Include(~bad_ab, [a, b]);
                            end if;
                        end for;
                    else
                        Include(~bad_cfs, cf);
                    end if;
                end if;
            end if;
        end if;
    end for;
    if # bad_cfs gt 0 then
        return false, [*bad_cfs, bad_ab*];
    end if;
    if # bad_ab eq 0 then
        return true, _;
    end if;
    bad_zwB := {};
    for ab in bad_ab do
        a := ab[1];
        b := ab[2];
        czsq := (a^2 + 4*b^2) * (a^2 + 2*a*b - 4*b^2)^2;
        cw := a^6 + 4*a^5*b + 64*a*b^5 - 64*b^6;
        c := Gcd(czsq, cw);
        zsq := Integers() ! (czsq / c);
        w := Integers() ! (cw / c);
        test, z := IsSquare(zsq);
        if test then
            DB4 := w^2 - z^4;
            if DB4 mod D eq 0 then
                B4 := Integers() ! (DB4 / D);
                test, B2 := IsSquare(B4);
                if test then
                    test, B := IsSquare(B2);
                    if test then
                        Include(~bad_zwB, [z, w, B]);
                    end if;
                end if;
            end if;
        end if;
    end for;
    ED := EllipticCurve([0, 0, 0, D, 0]);
    bad_points := {ED ! [zwB[1]^2*zwB[3], zwB[1]*zwB[2], zwB[3]^3] : zwB in bad_zwB};
    if bad_points eq { ED ! [0, 1, 0] } then
        return true, bad_points;
    else
        return false, bad_points;
    end if;
end function;

// We apply the function reducible_at_5 to various values of D
irreducible_at_5(-2);
// true { (0 : 1 : 0) }
irreducible_at_5(3);
// false { (1 : -2 : 1), (1 : 2 : 1), (0 : 1 : 0), (121/9 : -1342/27 : 1), (121/9 : 1342/27 : 1) }
irreducible_at_5(-17);
// true { (0 : 1 : 0) }
irreducible_at_5(125);
// true { (0 : 1 : 0) }

// Irreducible for l = 7
X014 := SmallModularCurve(14);
Genus(X014); // 1
j14 := jInvariant(X014, 14);
# Generators(X014); // 1
Q14 := Generators(X014)[1];
Order(Q14); // 6
{Evaluate(j14, d*Q14) : d in [1..6]}; // Infinity, -3375, 16581375
Factorization(Denominator(jE)); // t = 1, -1 corresponding to Infinity
Factorization(Numerator(jE) + 3375*Denominator(jE)); // t = -65/63 corresponding to -3375
Factorization(Numerator(jE) - 16581375*Denominator(jE)); // t = 65/63 corresponding to 16581375

// Irreducible for l = 11
X011 := SmallModularCurve(11);
Genus(X011); // 1
j11 := jInvariant(X011, 11);
# Generators(X011); // 1
Q11 := Generators(X011)[1];
Order(Q11); // 5
{Evaluate(j11, d*Q11) : d in [1..5]}; // -24729001, -32768, Infinity, -121
Factorization(Numerator(jE) + 24729001*Denominator(jE)); // no t corresponding to -24729001
Factorization(Numerator(jE) + 32768*Denominator(jE)); // no t corresponding to -32768
Factorization(Numerator(jE)); // t = -5/3 corresponding to Infinity
Factorization(Numerator(jE) + 121*Denominator(jE)); // no t corresponding to -121
