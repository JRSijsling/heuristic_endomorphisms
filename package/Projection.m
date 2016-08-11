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
gX := (Degree(fX) - 1) div 2;
KX<xX, yX> := FunctionField(X);
PX0 := X ! [ Coefficient(c, 0) : c in PX ];
duX := Differential(KX ! MyUniformizer(fX, PX0));
dxX := Differential(xX);
BX := [ dxX / yX , xX * dxX / yX ];
BXP := [ RP ! Evaluate(Q ! (b / duX), PX) : b in BX ];
pX := &+[ A[i] * BXP[i] : i in [1..#A] ];

fE := DefiningEquations(E)[1];
gE := (Degree(fE) - 1) div 2;
KE<xE, yE> := FunctionField(E);
PE0 := E ! [ Coefficient(c, 0) : c in PE ];
duE := Differential(KE ! MyUniformizer(fE, PE0));
dxE := Differential(xE);
BE := [ dxE / yE ];
BEP := [ RP !  Evaluate(Q ! (b / duE), PE) : b in BE ];
pE := BEP[1];

n := 1;
RP<t> := PowerSeriesRing(F, 2);
// Quotient always well-defined:
// (need absolute, not relative precision)
pE0 := Coefficient(pE, 0);
u := Coefficient(pX, 0)/pE0 * t + O(t^2);

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


function AlgebraizePoint(imPX, PX, X, deg);

prec := Precision(Parent(imPX[1])) - 1;
A2<x,y> := Ambient(X);
numdens := [ x^i : i in [0..(deg div 2)] ] cat [ x^i*y : i in [0..((deg - 3) div 2)] ];
A := Matrix([ [ Coefficient(Evaluate(numden, PX), i) : i in [0..(prec - 1)] ] : numden in numdens ]);
imPX_alg := [ ];
for c in imPX do
    B := Matrix([ [ Coefficient(Evaluate(numden, PX) * c, i) : i in [0..(prec - 1)] ] : numden in numdens ]);
    M := VerticalJoin(A, B);
    K := Kernel(M);
    w := K.1;
    num := [ w[i] : i in [1..#numdens] ];
    den := [ -w[i] : i in [(#numdens + 1)..(2*#numdens)] ];
    num := &+[ num[i] * numdens[i] : i in [1..#numdens] ];
    den := &+[ den[i] * numdens[i] : i in [1..#numdens] ];
    Append(~imPX_alg, num/den);
end for;
K<x,y> := FunctionField(X);
return [ K ! c : c in imPX_alg ];

end function;


function ProjectionToEllipticFactorG2(pX, pE, A, deg : margin := 2^4);
// Geared to a quite specific setting... for now

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

// TODO: Prec needed? Theory.
n := 4*deg + 2^4;
PX := DevelopInUnif(fX, PX0, n);
PE := DevelopInUnif(fE, PE0, n);
imPX := NthApproxsOld(X, PX, E, PE, A, n);
return AlgebraizePoint(imPX, PX, X, 2*deg);

end function;

// TODO: Optimize, good expressions, test
