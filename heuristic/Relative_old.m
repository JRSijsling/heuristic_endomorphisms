/***
 *  Relative number field functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


forward EmbedAtInfinitePlace;


intrinsic EmbedAsComplexPolynomials(fs::SeqEnum, prec::RngIntElt) -> SeqEnum
{Embeds a list of polynomials fs as complex polynomials to precision prec.}

R := Parent(fs[1]); d := #GeneratorsSequence(R);
F := BaseRing(R); iota := InfinitePlaces(F)[1];
CC<I> := ComplexField(prec); RCC := PolynomialRing(CC, d);

fsCC := [ EmbedAtInfinitePlace(f, iota, R, RCC) : f in fs ];
if d eq 1 then
    RCC := PolynomialRing(CC); fsCC := [ RCC ! fCC : fCC in fsCC ];
end if;

return fsCC, iota;

end intrinsic;


function EmbedAtInfinitePlace(f, iota, R, RCC)

if IsZero(f) then
    return RCC ! 0;
else
    prec := Precision(BaseRing(RCC));
    mons := Monomials(f);
    return &+[ Evaluate(MonomialCoefficient(f, mon), iota : Precision := prec) * Monomial(RCC, Exponents(mon)) : mon in mons ];
end if;

end function;


function CompareFields(K1, K2);
// Input:   Two subfields, fields, or polynomials.
// Output:  A comparison function: field with smaller degrees are smaller.

if Degree(K1) lt Degree(K2) then
    return -1;
elif Degree(K1) eq Degree(K2) then
    return 0;
else
    return 1;
end if;

end function;


intrinsic ClearFieldDenominator(K::FldNum) -> FldNum
{Simplifies the defining polynomial of a field to an integral version.}

F := BaseRing(K); d := Degree(K);
if d eq 1 then
    return F;
end if;
r := K.1;
coeffs := Coefficients(MinimalPolynomial(r, Rationals()));
dens := Reverse([ Denominator(coeff) : coeff in coeffs ]);
primes := &join[ Set([ tup[1] : tup in Factorization(den) ]) : den in dens | den ne 0 ];
if #primes eq 0 then
    common_den := 1;
else
    common_den := &*[ p^Maximum([ Ceiling(Valuation(dens[k], p)/k) : k in [1..#dens] | dens[k] ne 0 ]) : p in primes ];
end if;
f := MinimalPolynomial(common_den*r);
K := NumberField(f);

return K;

end intrinsic;


intrinsic ExtendRelativeSplittingField(K::FldNum, f::RngUPolElt) -> FldNum
{Extension step for relative splitting fields.}

print "Extending relative spl for", K, f;
F := BaseField(K);
if Degree(F) eq 1 then
    return SplittingField(f);
end if;
while true do
    factors := [ tup[1] : tup in Factorization(f, K) | Degree(tup[1]) gt 1 ];
    if #factors eq 0 then
        return K;
    end if;
    print "Current F:", F, BaseRing(F);
    print "Current K:", K, BaseRing(K);
    print IsSubfield(F, NumberField(factors[1]));
    print MinimalPolynomial(K.1, F);
    print NumberField(factors[1]);
    K := RelativeField(F, NumberField(factors[1]));
    print "After RelativeField:", K;
end while;

end intrinsic;


intrinsic RelativeSplittingField(fs::SeqEnum) -> FldNum
{Determines a relative splitting field of the polynomials in fs.}

F := BaseRing(fs[1]);
R<x> := PolynomialRing(F);
K := NumberField(x : DoLinearExtension := true);
fs := Reverse(Sort(fs, CompareFields));
for f in fs do
    if not HasRoot(f, K) then
        for tup in Factorization(f, K) do
            K := ExtendRelativeSplittingField(K, tup[1]);
            return K;
            K := ClearFieldDenominator(K);
        end for;
    end if;
end for;
return K;

end intrinsic;


intrinsic RelativeSplittingField(f::RngUPolElt) -> FldNum
{Determines a relative splitting field of the polynomials f.}

return RelativeSplittingField([ f ]);

end intrinsic;


// recognition, overfield, lattice, data
