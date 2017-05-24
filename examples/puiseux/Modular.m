AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);

F := Rationals();
R<x> := PolynomialRing(F);
f := x^8 - 12*x^7 + 50*x^6 - 108*x^5 + 131*x^4 - 76*x^3 - 10*x^2 + 44*x - 19;
X := HyperellipticCurve(f);
P0 := X ! [1, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[-1, -2, -4],
[ 2,  3,  4],
[-1, -1, -1]
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
