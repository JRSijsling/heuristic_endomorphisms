/***
 *  Relative number field functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function EmbedAsComplexPolynomials(f, h : prec := prec0);

CC<I> := ComplexField(prec);
RCC<xCC> := PolynomialRing(CC);
R<x> := Parent(f);
F := BaseRing(R);
infs := InfinitePlaces(F);
fCC := &+[ Evaluate(Coefficient(f, i), infs[1] : Precision := prec) * xCC^i : i in [0..Degree(f)] ];
if IsZero(h) then
    hCC := RCC ! 0;
else
    hCC := &+[ Evaluate(Coefficient(h, i), infs[1] : Precision := prec) * xCC^i : i in [0..Degree(h)] ];
end if;
iota := Evaluate(F.1, infs[1] : Precision := prec);

return fCC, hCC, iota;

end function;

