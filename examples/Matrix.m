R<x> := PolynomialRing(Rationals());
K<r> := NumberField(x^2 + 2);
M := MatrixRing(K, 3);

AK := AssociativeAlgebra(M);
OOK := MaximalOrder(AK);
//print Trivialize(AK, Basis(OOK));

AQQ := RestrictionOfScalars(AK);
OOQQ := MaximalOrder(AQQ);
//print Trivialize(AQQ, Basis(OOQQ));

exit;
