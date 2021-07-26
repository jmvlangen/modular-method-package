/* Conductor exponent at 2 for E_1, case c odd */
mytime := Cputime();

p := 2;
M := p^6;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 0, 0, -3*2*(s^4 + 2*t*s^3 + 2*t^3*s + t^4), -2*3*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)]);
            LI := LocalInformation(E, p);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.370 s

/* Conductor exponent at 3 for E_1, case c odd */
mytime := Cputime();

p := 3;
M := p^6;
listp := [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 0, 0, -3*2*(s^4 + 2*t*s^3 + 2*t^3*s + t^4), -2*3*(s - t)*(s + t)*(s^4 + 2*s^3*t + 6*s^2*t^2 + 2*s*t^3 + t^4)]);
            LI := LocalInformation(E, p);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 124.490 s

/* Conductor exponent at prime above 2 for E_2, case c odd */
mytime := Cputime();

p := 2;
K<sqrt3> := QuadraticField(3);
O := MaximalOrder(K);
L := Factorization(2*O);
P := L[1][1];
M := p^3;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 2*(sqrt3 - 1)*(s - t), 0, (2 - sqrt3)*((s - t)^2 - 2*sqrt3*s*t), 0]);
            LI := LocalInformation(E, P);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),P),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.020 s

/* Conductor exponent at 2 for E_1, case c even, case (3.2) */
mytime := Cputime();

p := 2;
M := p^6;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 0, 0, -3*(3*s^4 + 6*t^2*s^2 - t^4), -2*(6*s*t*(3*s^4 + t^4))]);
            LI := LocalInformation(E, p);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.360 s

/* Conductor exponent at 2 for E_1, case c even, case (3.3) */
mytime := Cputime();

p := 2;
M := p^6;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 0, 0, -12*(-3*s^4 + 6*t^2*s^2 + t^4), -16*(6*s*t*(3*s^4 + t^4))]);
            LI := LocalInformation(E, p);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.390 s

/* Conductor exponent at 3 for E_1, case c even, case (3.2) */
mytime := Cputime();

p := 3;
M := p^6;
listp := [[0, 1], [0, 2], [1, 1], [1, 2], [2, 1], [2, 2]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 0, 0, -3*(3*s^4 + 6*t^2*s^2 - t^4), -2*(6*s*t*(3*s^4 + t^4))]);
            LI := LocalInformation(E, p);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 109.370 s

/* Conductor exponent at 3 for E_1, case c even, case (3.3) */
mytime := Cputime();

p := 3;
M := p^6;
listp := [[0, 1], [0, 2], [1, 1], [1, 2], [2, 1], [2, 2]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 0, 0, -12*(-3*s^4 + 6*t^2*s^2 + t^4), -16*(6*s*t*(3*s^4 + t^4))]);
            LI := LocalInformation(E, p);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 124.220 s

/* Conductor exponent at prime above 2 for E_2, case c even, positive sign */
mytime := Cputime();

p := 2;
K<sqrt3> := QuadraticField(3);
O := MaximalOrder(K);
L := Factorization(2*O);
P := L[1][1];
M := p^3;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 4*(sqrt3 - 1)*t, 0, (sqrt3 - 1)^2*(sqrt3*s^2 + (-2 + sqrt3)*t^2), 0]);
            LI := LocalInformation(E, P);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),P),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.030 s

/* Conductor exponent at prime above 2 for E_2, case c even, negative sign */
mytime := Cputime();

p := 2;
K<sqrt3> := QuadraticField(3);
O := MaximalOrder(K);
L := Factorization(2*O);
P := L[1][1];
M := p^3;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 4*(sqrt3 - 1)*t, 0, (sqrt3 - 1)^2*(sqrt3*s^2 + (-2 - sqrt3)*t^2), 0]);
            LI := LocalInformation(E, P);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),P),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.030 s

/* Conductor exponent at prime above 2 for E_3, case c even, positive sign */
mytime := Cputime();

p := 2;
K<sqrt3> := QuadraticField(3);
O := MaximalOrder(K);
L := Factorization(2*O);
P := L[1][1];
M := p^3;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 + 3)*s^2), 0]);
            LI := LocalInformation(E, P);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),P),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.030 s

/* Conductor exponent at prime above 2 for E_3, case c even, negative sign */
mytime := Cputime();

p := 2;
K<sqrt3> := QuadraticField(3);
O := MaximalOrder(K);
L := Factorization(2*O);
P := L[1][1];
M := p^3;
listp := [[0, 1], [1, 0]];

for s := 0 to M - 1 do
    for t := 0 to M - 1 do
        if ([s mod p, t mod p] in listp) then
            E := EllipticCurve([0, 12*(sqrt3 - 1)*s, 0, 3*sqrt3*(sqrt3 - 1)^2*(t^2 + (2*sqrt3 - 3)*s^2), 0]);
            LI := LocalInformation(E, P);
            printf "%3o %3o %3o %3o %3o %3o\n",s,t,LI[3],LI[5],Valuation(Discriminant(E),P),LI[2];
        end if;
    end for;
end for;

mytime := Cputime() - mytime;
printf "Time elapsed: %3o",mytime; // 0.040 s
