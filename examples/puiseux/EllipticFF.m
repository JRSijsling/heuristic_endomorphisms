AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

F := FiniteField(17, 2);
R<t> := PolynomialRing(F);
q := t^2 + t + 3;
r := Roots(q)[1][1];

R<x> := PolynomialRing(F);
f := x^3 - x^2 - 7*x + 10;
h := 1;
p := 4*f + h^2;
X := HyperellipticCurve(p);
P0 := X ! [2, 1];
//P0 := X ! [1, 0, 0];

M := Matrix(F, [[2]]);
//M := Matrix(F, [[r]]);

print "Field:";
print F;
print "Curve:";
print X;

print "Tangent representation:";
print M;
print "Minimal polynomial:";
print MinimalPolynomial(M);

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

exit;
