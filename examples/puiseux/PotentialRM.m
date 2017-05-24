AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - t - 1);
R<x> := PolynomialRing(F);
f := x^6 + 2*x^3 - x;
h := x^3 + 1;
p := 4*f + h^2;
X := HyperellipticCurve(p);
print X;
P0 := X ! [0, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[   -r,     0],
[    0, r - 1]
]);

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
time test, fs := CantorMorphismFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print fs;

exit;
