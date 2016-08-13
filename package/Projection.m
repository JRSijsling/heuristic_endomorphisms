function MyUniformizer(f, P0);

R<x, y> := Parent(f);
if Evaluate(Derivative(f, y), P0) ne 0 then
    return x;
else
    return y;
end if;

end function;


function MyEvaluate(f, P, tj);

R<x, y> := Parent(f);
xy := [x, y];
F := BaseRing(R);
RP<t> := Parent(tj);
S := PolynomialRing(RP);
P0 := [ Coefficient(c, 0) : c in P ];
u := MyUniformizer(f, P0);
if u eq x then
    i := 1;
else
    i := 2;
end if;
subst := [S!0, S!0];
subst[i] := tj;
subst[(i mod 2) + 1] := S.1;
h := hom<R -> S | subst>;
g := h(f);
for rt in [ tup[1] : tup in Roots(g) ] do
    if Coefficient(rt, 0) eq Coefficient(P[(i mod 2) + i], 0) then
        dev := [ RP!0, RP!0 ];
        dev[i] := tj;
        dev[(i mod 2) + 1] := rt;
        return dev;
    end if;
end for;

end function;


function MySum(L);

if #L eq 0 then
    return 0;
end if;
return &+L;

end function;


function DevelopInUnif(f, P0, n);

R<x, y> := Parent(f);
xy := [x, y];
F := BaseRing(R);
RP<t> := PowerSeriesRing(F, n);
u := MyUniformizer(f, P0);
S := PolynomialRing(RP);
if u eq x then
    i := 1;
else
    i := 2;
end if;
subst := [S!0, S!0];
subst[i] := t + P0[i];
subst[(i mod 2) + 1] := S.1;
h := hom<R -> S | subst>;
g := h(f);
for rt in [ tup[1] : tup in Roots(g) ] do
    if Coefficient(rt, 0) eq P0[(i mod 2) + i] then
        dev := [ RP!0, RP!0 ];
        dev[i] := t + P0[i];
        dev[(i mod 2) + 1] := rt;
        return dev;
    end if;
end for;

end function;


function NthApproxsOld(X, PX, E, PE, A, N);

A2 := Ambient(X);
Q<x,y> := FieldOfFractions(CoordinateRing(A2));
F := BaseRing(X);
RP<t> := PowerSeriesRing(F);

fX := DefiningEquations(X)[1];
KX<xX, yX> := FunctionField(X);
PX0 := X ! [ Coefficient(c, 0) : c in PX ];
duX := Differential(KX ! MyUniformizer(fX, PX0));
dxX := Differential(xX);
BX := [ dxX / yX , xX * dxX / yX ];
BXP := [ RP ! Evaluate(Q ! (b / duX), PX) : b in BX ];
pX := &+[ A[i] * BXP[i] : i in [1..#A] ];

fE := DefiningEquations(E)[1];
KE<xE, yE> := FunctionField(E);
PE0 := E ! [ Coefficient(c, 0) : c in PE ];
duE := Differential(KE ! MyUniformizer(fE, PE0));
dxE := Differential(xE);
BE := [ dxE / yE ];
BEP := [ RP !  Evaluate(Q ! (b / duE), PE) : b in BE ];
pE := BEP[1];

n := 1;
RP<t> := PowerSeriesRing(F, 2);
pE0 := Coefficient(pE, 0);
u := (Coefficient(pX, 0)/pE0)*t + O(t^2);

while n lt N do
    n := n + 1;
    RP<t> := PowerSeriesRing(F, n + 1);
    u := &+[ Coefficient(u, j) * t^j : j in [1..n-1] ] + O(t^(n + 1));
    ev_lhs := Evaluate(RP ! pE, u) * Derivative(u);
    ev_rhs := RP ! pX;
    un := Coefficient(ev_rhs - ev_lhs, n - 1) / (pE0 * n);
    u +:= un * t^n;
end while;

return MyEvaluate(fE, PE, u);

end function;


function LinTermMult(lt1, lt2);

n := #Rows(Transpose(lt1));
return Matrix([ [ &+[ lt1[1, j] * lt2[1, i + 1 - j] : j in [1..i] ] : i in [1..n] ] ]);

end function;


function EvaluateOldGuess(fi, tj_old, n, m);

F := BaseRing(Parent(tj_old));
RP<t> := PowerSeriesRing(F, n + m + 1);
tj := &+[ Coefficient(tj_old, i) * t^i : i in [1..n] ] + O(t^(n + m + 1));

power0 := RP ! 1;
power1 := tj;
power_term0 := Matrix([ [ F ! 0 : k in [1..m] ] ]);
power_term1 := Matrix([ [ F ! 1 ] cat [ F ! 0 : k in [2..m] ] ]);
deriv := Matrix([ [ Coefficient(tj, k) : k in [1..m] ] ]);

powers := [ power0, power1 ];
power_terms := [ power_term0, power_term1 ];
power := power1;
power_term := power_term1;
for k in [2..(n + m - 1)] do
    Append(~powers, power1 * powers[#powers]);
    if k le (m - 1) then
        power_term := LinTermMult(deriv, power_term);
        Append(~power_terms, k * power_term);
    end if;
end for;

evs0 := [ ];
lin_terms := [ ];
ev0 := &+[ Coefficient(fi, k) * powers[k + 1] : k in [0..(n + m - 1)] ] + O(t^(n + m));
lin_term := Coefficient(fi, 0) * power_terms[1];
for k:=2 to m do
    c := Coefficient(fi, k - 1);
    power_term := power_terms[k];
    for l:=k to m do
        lin_term[1, l] +:= c * power_term[1, l + 1 - k];
    end for;
end for;

return ev0, lin_term;

end function;


function NthApproxs(X, PX, E, PE, A, N);

A2 := Ambient(X);
Q<x,y> := FieldOfFractions(CoordinateRing(A2));
F := BaseRing(X);
RP<t> := PowerSeriesRing(F);

fX := DefiningEquations(X)[1];
KX<xX, yX> := FunctionField(X);
PX0 := X ! [ Coefficient(c, 0) : c in PX ];
duX := Differential(KX ! MyUniformizer(fX, PX0));
dxX := Differential(xX);
BX := [ dxX / yX , xX * dxX / yX ];
BXP := [ RP ! Evaluate(Q ! (b / duX), PX) : b in BX ];
pX := &+[ A[i] * BXP[i] : i in [1..#A] ];

fE := DefiningEquations(E)[1];
KE<xE, yE> := FunctionField(E);
PE0 := E ! [ Coefficient(c, 0) : c in PE ];
duE := Differential(KE ! MyUniformizer(fE, PE0));
dxE := Differential(xE);
BE := [ dxE / yE ];
BEP := [ RP ! Evaluate(Q ! (b / duE), PE) : b in BE ];
pE := BEP[1];

n := 1;
RP<t> := PowerSeriesRing(F, 2);
pE0 := Coefficient(pE, 0);
u := (Coefficient(pX, 0)/pE0)*t + O(t^2);
while n lt N do
    m := Minimum(n, N - n);
    RP<t> := PowerSeriesRing(F, n + m + 1);
    ev0, lin_term := EvaluateOldGuess(pE, u, n, m);
    u := &+[ Coefficient(u, j) * t^j : j in [1..n] ] + O(t^(n + m + 1));
    ev_lhs := ev0 * Derivative(u);
    ev_rhs := RP ! pX;
    diff1 := ev_rhs - ev_lhs;
    for k in [1..m] do
        diff2 := Coefficient(diff1, n + k - 1);
        for l in [1..k] do
            diff2 -:= (n + l) * Coefficient(u, n + l) * Coefficient(ev0, k - l);
            gamma := MySum([ lin_term[1, l + 1 - ind] * Coefficient(u, n + ind) : ind in [1..Minimum(l, k - 1)] ]);
            diff2 -:= (k + 1 - l) * Coefficient(u, k + 1 - l) * gamma;
        end for;
        un := diff2 / (pE0 * (n + k));
        u +:= un * t^(n + k);
    end for;
    n *:= 2;
end while;

return MyEvaluate(fE, PE, u);

end function;


function AlgebraizeTangentRepresentation(X, fX, PX0, E, fE, PE0, A, degf0 : margin := 2^4);

for degf in [degf0..degf0^2] do
    print degf;
    deg := 3*degf + 3;
    n := 2*deg + margin;

    repeat
        PX := DevelopInUnif(fX, PX0, n);
        PE := DevelopInUnif(fE, PE0, n);
        imPX := NthApproxs(X, PX, E, PE, A, n);

        prec := Precision(Parent(imPX[1])) - 1;
        A2<x,y> := Ambient(X);
        numdens := [ x^i : i in [0..(deg div 2)] ] cat [ x^i*y : i in [0..((deg - 5) div 2)] ];
        A1 := Matrix([ [ Coefficient(Evaluate(numden, PX), i) : i in [0..(prec - 1)] ] : numden in numdens ]);
        n +:= 1;
    until Rank(A1) eq #numdens;

    imPX_alg := [ ];
    for c in imPX do
        A2 := Matrix([ [ Coefficient(Evaluate(numden, PX) * c, i) : i in [0..(prec - 1)] ] : numden in numdens ]);
        M := VerticalJoin(A1, A2);
        K := Kernel(M);
        if Dimension(K) ne 0 then
            w := K.1;
            num := [ w[i] : i in [1..#numdens] ];
            den := [ -w[i] : i in [(#numdens + 1)..(2*#numdens)] ];
            num := &+[ num[i] * numdens[i] : i in [1..#numdens] ];
            den := &+[ den[i] * numdens[i] : i in [1..#numdens] ];
            Append(~imPX_alg, num/den);
        end if;
    end for;

    if #imPX_alg eq 2 then
        KX<x,y> := FunctionField(X);
        eqs := [ KX ! c : c in imPX_alg ];
        if Evaluate(fE, eqs) eq 0 then
            return true, eqs;
        end if;
    end if;
end for;

return false, 0;

end function;


function ProjectionToEllipticFactorG2(pX, pE, A, deg0 : margin := 2^4);
// Geared to a quite specific setting for now

S<t> := Parent(pX);
F := SplittingField((t^2 - Evaluate(pX, 0))*(t^2 - Evaluate(pE, 0)));
S<t> := PolynomialRing(F);
pX := S ! pX;
pE := S ! pE;
A := [ F ! c : c in A ];

A2 := AffineSpace(F, 2);
R<x,y> := CoordinateRing(A2);
h := hom<S -> R | x>;
fX := y^2 - h(pX);
fE := y^2 - h(pE);
X := Curve(A2, fX);
E := Curve(A2, fE);
KX := FunctionField(X);
KE := FunctionField(E);
PX0 := X ! [ 0, Roots(t^2 - Evaluate(pX, 0))[1][1] ];
PE0 := E ! [ 0, Roots(t^2 - Evaluate(pE, 0))[1][1] ];

return AlgebraizeTangentRepresentation(X, fX, PX0, E, fE, PE0, A, deg0);

end function;


function MyTest();

n := 4;
m := 4;
RP<x> := PowerSeriesRing(Rationals(), 2*n + 1);
fi := &+[ l*x^(l-1) : l in [1..(2*n + 1)] ];
K<a,b,c,d> := RationalFunctionField(Rationals(), n);
RP<t> := PowerSeriesRing(K, 2*n);
tj := 111*t^1 + 17*t^2 - 11*t^3 + 31*t^4 + a*t^5 + b*t^6 + c*t^7 + d*t^8;
ev0, M := EvaluateOldGuess(fi, tj, n, m);

return [* [ tj^k : k in [1..2*n] ], Evaluate(fi, tj), ev0, M *];

end function;
