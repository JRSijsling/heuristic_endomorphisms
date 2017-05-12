K<mu> := RationalFunctionField(Rationals());
R<x> := PolynomialRing(K);

f := -x + mu^2;
g := x^2 + mu^2*x + mu^4;
h := mu^6 + 1;

f := R ! f;
g := R ! g;
h := R ! h;

print h;
print f*g;

f2 := Coefficient(f, 2); f1 := Coefficient(f, 1); f0 := Coefficient(f, 0);
g2 := Coefficient(g, 2); g1 := Coefficient(g, 1); g0 := Coefficient(g, 0);
h2 := Coefficient(h, 2); h1 := Coefficient(h, 1); h0 := Coefficient(h, 0);

A := Matrix([ [ f2, f1, f0 ], [ h2, h1, h0 ], [ g2, g1, g0 ] ]);
print A;
B := A^(-1);
print B;

a1 := B[1,1]; b1 := B[1,2]; c1 := B[1,3];
a2 := B[2,1]; b2 := B[2,2]; c2 := B[2,3];
a3 := B[3,1]; b3 := B[3,2]; c3 := B[3,3];

a := a1 + 2*a2*x + a3*x^2;
b := b1 + 2*b2*x + b3*x^2;
c := c1 + 2*c2*x + c3*x^2;

f := 2*b*(b^2 - a*c);
print f;
print ((mu^18 + 3*mu^12 + 3*mu^6 + 1)/2)*f;
X := HyperellipticCurve(f);

exit;

