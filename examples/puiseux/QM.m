AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 + 3);
R<x> := PolynomialRing(F);
f := 6*x^5 + 9*x^4 - x^3 - 3*x^2;
h := 1;
p := 4*f + h^2;
X := HyperellipticCurve(p);
print X;
P0 := X ! [0, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[   -r,     r],
[  2*r,     r]
]);
M := Transpose(M);

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
