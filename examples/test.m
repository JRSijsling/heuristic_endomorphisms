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

exit;
