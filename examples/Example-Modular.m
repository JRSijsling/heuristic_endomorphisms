AttachSpec("../spec");
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

time fs := CantorMorphismFromMatrixSplit(X, P0, M : LowerBound := 1);
time D := DivisorFromMatrixSplit(X, P0, M : LowerBound := 1);

DEs := DefiningEquations(D);
DEs := [ 2*30697399*DEs[1], 2*30697399*DEs[2] ];
Lat := Lattice(Matrix([ Coefficients(DE) : DE in DEs ]));
print LLL(Lat);

exit;
