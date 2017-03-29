/***
 *  Canonizing a matrix
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function CanonizeMatrices(grep, frep, idems);

R<x> := PolynomialRing(Rationals());
f := R ! frep;
// FIXME: Magma greatly loses here since the rationals do not admit a subfield constructor
if frep eq [0, 1] then
    L := Rationals();
    K := Rationals();
else
    L := NumberField(f);
    K := sub< L | L ! grep >;
end if;

return [ Matrix(K, [ [ K ! Eltseq(c) : c in Eltseq(row) ] : row in Rows(idem) ]) : idem in idems ];

end function;
