/***
 *  Relative number field functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function EmbedAsComplexPolynomials(fs : prec := prec0);

R := Parent(fs[1]); F := BaseRing(R); infs := InfinitePlaces(F);
CC<I> := ComplexField(prec);
iota := Evaluate(F.1, infs[1] : Precision := prec);

gens := GeneratorsSequence(R); d := #gens;
RCC := PolynomialRing(CC, d); gensCC := GeneratorsSequence(RCC);
fsCC := [ EmbedAsComplexPolynomial(f, infs[1], R, RCC) : f in fs ];
// FIXME: Coerce to univariate because Magma is idiotic
if d eq 1 then
    RCC := PolynomialRing(CC);
end if;

return [ RCC ! f : f in fs ], iota;

end function;


function EmbedAsComplexPolynomial(f, inf, R, RCC);

if IsZero(f) then
    return RCC ! 0;
else
    prec := Precision(BaseRing(RCC));
    mons := Monomials(f);
    return &+[ Evaluate(MonomialCoefficient(f, mon), inf : Precision := prec) * Monomial(RCC, Exponents(mon)) : mon in mons ];
end if;

end function;

