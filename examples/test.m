AttachSpec("../spec");

R<x> := PolynomialRing(Rationals());
X := HyperellipticCurve(x^5 + x + 1);
print NonWeierstrassBasePoint(X, Rationals());
print NonWeierstrassBasePoint(X, NumberField(x^2 - 5));

P2 := ProjectiveSpace(Rationals(), 2);
P2<x,y,z> := ProjectiveSpace(Rationals(), 2);
f := x^3*y + y^3*z + z^3*x;
X := Curve(Scheme(P2, f));
P := NonWeierstrassBasePoint(X, Rationals());
print Parent(P[1]);

R<x> := PolynomialRing(Rationals());
F<r> := NumberFieldExtra(x^2 - 2);
R<x> := PolynomialRing(F);
K<s> := NumberFieldExtra(x^2 - (1 + r));
print K;

prec := 500; CC<I> := ComplexFieldExtra(prec);
aCC := CC ! Evaluate(s, K`iota : Precision := prec);
print RelativeMinimalPolynomial(aCC, F);
print AlgebraizeElementInRelativeField(aCC, K);

print Eltseq(s);
print Eltseq(Eltseq(s)[1]);

R<x> := PolynomialRing(Rationals());
pols :=
[
x^2 + 1, x^2 + 4, x^2 + 1, x^2 + 1, x^2 + 2, x^2 + 1/2, x^2 - 2, x^2 - 2, x^2 - 1/2, x^2 - 2, x - 1, x, x, x - 1, x, x
];

print pols;
print RelativeSplittingFieldExtra(pols);

exit;
