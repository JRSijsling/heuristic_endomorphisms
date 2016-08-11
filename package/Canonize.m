function CanonizeMatrix(grep, frep, idem);

R<x> := PolynomialRing(Rationals());
f := R ! frep;
// TODO: Magma greatly loses here since the rationals do not admit a subfield constructor
if frep eq [0, 1] then
    L := Rationals();
    K := Rationals();
else
    L := NumberField(f);
    K := sub< L | L ! grep >;
end if;

return Matrix(K, [ [ K ! c : c in Eltseq(row) ] : row in Rows(idem) ]);

end function;
