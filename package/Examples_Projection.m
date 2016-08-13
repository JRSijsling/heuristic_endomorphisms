load "Projection.m";

R<x> := PolynomialRing(Rationals());

fX := x^6 + x^4 + x^2 + 1;
fE := x^3 + x^2 + x + 1;
A := [0, 2];
deg := 2;

fX := x^6 + x^4 + x^2 + 1;
fE := x^4 + x^3 + x^2 + x;
A := [0, 2];
deg := 2;

fX := x^5 + x^3 + x;
fE := x*(x - 2)*(x + 1);
A := [-1, 1];
deg := 2;

fX := x^5 + x^3 + x;
fE := (x + 1)*(x - 1)*(x + 2);
A := [-1, 1];
deg := 2;

print ProjectionToEllipticFactorG2(fX, fE, A, deg : margin := 100);

exit;
