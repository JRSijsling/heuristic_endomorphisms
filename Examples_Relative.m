AttachSpec("spec");

F := Rationals();
R<x> := PolynomialRing(F);
f := x^3 + x + 1;

K := NumberField(f);
print K;
print BaseRing(K);

L := RelativeSplittingField([f, x]);
print L;
print BaseRing(L);
print Roots(f, L);

F<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(F);
f := x^3 + x - r;

K := NumberField(f);
print K;
print BaseRing(K);

L := RelativeSplittingField([f, x]);
print L;
print BaseRing(L);
print Roots(f, L);

print InfinitePlaces(F);
print InfinitePlaces(L);

infF := InfinitePlaces(F)[1];
for infL in InfinitePlaces(L) do
    print Extends(infL, infF);
end for;

exit;
