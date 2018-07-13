for p in PrimesInInterval(5,100) do         
F := Factorization(p*OK);
if #F eq 1 then
P:=F[1,1];
# Reduction(E,P);
for f in nfs0 do
Norm(eps(p) * p + 1 - Coefficient(f[1], p));
end for;
"*";
end if;
end for;

for p in PrimesInInterval(5,100) do         
F := Factorization(p*OK);
if #F eq 2 then
    P:=F[1,1];
    Q:=F[2,1];
#Reduction(E,P) * #Reduction(E,Q);
for f in nfs0 do
Norm(eps(p) * p + 1 - Coefficient(f[1], p));
end for;
"*";
end if;
end for;

for p in PrimesInInterval(5,1000) do         
F := Factorization(p*OK);
if #F eq 2 then
    P:=F[1,1];
    Q:=F[2,1];
#Reduction(E,P) * #Reduction(E,Q) eq Norm(eps(p) * p + 1 - Coefficient(nfs[3,1], p));
end if;
end for;


