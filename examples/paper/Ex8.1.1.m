AttachSpec("../../spec");
SetVerbose("EndoCheck", 0);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - t - 1);
R<x> := PolynomialRing(F);
p := 5*x^6 + 10*x^3 - 4*x + 1;
X := HyperellipticCurve(p);
P0 := X ! [0, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[ r, 0     ],
[ 0, 1 - r ]
]);
M := -M;

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print fs;

eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
S<x2,x1> := PolynomialRing(F, 2);
res := hom<R -> S | [0, 0, x2, x1]>;
I := ideal<S | res(eqs[#eqs])>;
G := GroebnerBasis(I);

A := AffineSpace(S);
D := Scheme(A, G);
Is := IrreducibleComponents(D);
print "Divisor on PP^1 x PP^1:";
print Is[1];

exit;
