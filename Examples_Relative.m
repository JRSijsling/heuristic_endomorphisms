AttachSpec("spec");

F := Rationals();
R<x> := PolynomialRing(F);
f := x^3 + x + 1;
g := x^2 - 5;

K := NumberFieldExtra(f);
SetInfinitePlace(K, InfinitePlaces(K)[2]);

K := NumberField(f);
L := RelativeSplittingFieldExtra([f, g]);
print L;
print L`iota;
//print Roots(f, L);

F<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(F);
f := x^3 + x - r;
g := x^2 - 5;
g := x^2 - 5 + r;

K := NumberField(f);
L := RelativeSplittingFieldExtra([f, g]);
print L;
print L`iota;
//print Roots(f, L);

SetInfinitePlace(L, InfinitePlaces(L)[7]);
print L`iota;
print F`iota;

exit;
