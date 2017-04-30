AttachSpec("../spec");

R<x> := PolynomialRing(Rationals());
num := 101*x; den := 1;
num := x + 3; den := 5*x + 2;
num := 2*x + 3; den := 4*x + 5;
f := R ! (den^6 * Evaluate(x^5 + x + 1, num/den));
X := HyperellipticCurve(f);
g2X := G2Invariants(X);

prec := 300;
eqsCC := EmbedCurveEquations(X, prec);
P := PeriodMatrix(eqsCC : HaveOldenburg := true);
K := BaseRing(X);

print X;
Y := FactorReconstructG2(P, K);
print Y;

exit;

print X;
print Y;
print g2X;
print g2Y;

exit;

R<x> := PolynomialRing(Rationals());
K<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(K);
X := HyperellipticCurve(x^5 + r*x + 1, x);
X := HyperellipticCurve(x^5 + 2*x + r, x);
P := NonWeierstrassBasePoint(X, K);
print P;
print MinimalPolynomial(P[2]);

exit;

K := Rationals();
R<x> := PolynomialRing(K);
X := HyperellipticCurve(x^3 + x + 1);
X := HyperellipticCurve(x^6 + 1);
X := HyperellipticCurve(2*x^6 + 1);
X := HyperellipticCurve(2*x^6 + 3*x + 7);
X := HyperellipticCurve(2*x^5 + 3*x + 7);
X := HyperellipticCurve((x^5 + x + 1)/12);
P := NonWeierstrassBasePoint(X, K);
print P;
print MinimalPolynomial(P[2]);

exit;

R<x> := PolynomialRing(Rationals());
K<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(K);
X := HyperellipticCurve(x^5 + r*x + 1, x);
P2<x,y,z> := ProjectiveSpace(K, 2);
X := Curve(Scheme(P2, x*y^3 + y*z^3 + r*z*x^3));
print X;
print FactorDescription(X);

exit;
