AttachSpec("spec");

CC := ComplexField(300);
R<x> := PolynomialRing(CC);

f := x^5 + x^4 + 2*x^3 + x^2 + x;
h := x^2 + x;
f := x^7 + x^6 + x^5 + x^3 + x^2 + x;
h := x^4 + x^2 + 1;

g := 4*f + h^2;
print PeriodMatrix(g);

exit;

