R<x> := PolynomialRing(Rationals());
F<r> := NumberField(x^3 + x + 1);

R<x> := PolynomialRing(F);
f1:= x^4 + x + 1;
K1<s1> := NumberField(f1);

R<x> := PolynomialRing(K1);
f2:= (x^4 + x + 1) div (x - s1);
K2<s2> := NumberField(f2);

R<x> := PolynomialRing(K2);
f3:= (x^4 + x + 1) div ((x - s1)*(x - s2));
K3<s3> := NumberField(f3);

print F;
print K1;
print K2;
print K3;

L := AbsoluteField(K3);
K := RelativeField(F, L);
print K;

exit;
