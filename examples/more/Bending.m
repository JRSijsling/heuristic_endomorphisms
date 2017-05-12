R<x> := PolynomialRing(Rationals());
F<u> := NumberField(x^2 + 6);

P2<x,y,z> := ProjectiveSpace(F, 2);
C := Conic(P2, x^2 + 3*y^2 + z^2);
P := C ! [1, u/3, 1];
phi := Parametrization(C, P);
P1 := Domain(phi);
C := Codomain(phi);

// Parameter:
t := 2 - u;
P := P1 ! [t, 1];
Q := phi(P);
print Q;

a,b,c := Explode(Eltseq(Q));
sqrt1 := (4/3)*(b/a);
sqrt2 := 4*(c/a);
j := sqrt1^2;
print j;

R<x> := PolynomialRing(F);
s := u*sqrt1;
t := 2*sqrt2^2;

f := (-4 + 3*s)*x^6 + 6*t*x^5 + 3*t*(28 + 9*s)*x^4 - 4*t^2*x^3 + 3*t^2*(28 - 9*s)*x^2 + 6*t^3*x - t^3*(4 + 3*s);
X := HyperellipticCurve(f);
I := G2Invariants(X);

print s;
print t;
print X;
print I;

exit;

