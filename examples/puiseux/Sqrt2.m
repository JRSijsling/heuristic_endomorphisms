AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

F := Rationals();
R<t> := PolynomialRing(F);
f := x^6 + 2*x^5 + 7*x^4 + 6*x^3 + 13*x^2 + 4*x + 8;
X := HyperellipticCurve(2*f);
P0 := X ! [0, 4, 1];

M := Matrix(F, [
[0, 1],
[2, 0]
]);
M := 1 + M;
//M := M^2;

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

exit;
