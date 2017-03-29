/***
 *  Elliptic curves from decompositions
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function DecompositionDegree(R);

return LCM([ Denominator(c) : c in Eltseq(R) ]);

end function;


function SmallIdempotent(Q, GensQ);
// Input:   A split quaternion algebra and a basis for an order.
// Output:  A nilpotent element that is small with respect to that order in the
//          sense that the numerator is small

test, M, f := IsMatrixRing(Q : Isomorphism := true);
invf := Inverse(f);

// TODO: What follows is obviously stupid (it does not even use GensQ) and
// should be improved.
return invf(M![1,0,0,0]);

end function;


function InducedEmbedding(fsubgen, fhom, frep : epscomp := epscomp0);
// Returns induced embedding on a subfield.

CC := Parent(fhom);
subhom := CC ! 0;
power := CC ! 1;
for i := 1 to #fsubgen do
    subhom := subhom + fsubgen[i] * power;
    power := power * fhom;
end for;

// Insert error test to be certain:
L := NumberField(Polynomial(frep));
s := L ! fsubgen;
p := MinimalPolynomial(s);
if not Abs(Evaluate(p, subhom)) lt epscomp then
    error Error("Failure in computation of induced embedding of subfield");
end if;

return subhom;

end function;


function LatticesFromIdempotents(idemsAn, P : epscomp := epscomp0, epsLLL :=
    epsLLL0, epsinv := epsinv0);
// Input:   Representations of idempotents terms of an endomorphism basis, that
//          basis in its analytic representation, and the corresponding period
//          matrix.
// Output:  Lattices corresponding to the elliptic curve quotients.

Ls := [ ]; col_numbers := [ ];
for idem in idemsAn do
    // Create analytic idempotent and project:
    PEllHuge := P * idem;
    col_number := 0;
    repeat
        col_number +:= 1;
        PEllBig := Transpose(Matrix([ [ row[col_number] : row in Rows(PEllHuge) ] ]));
    until not &and([ Abs(c) lt epscomp : c in Eltseq(PEllBig) ]);
    PEllBigSplit := SplitPeriodMatrix(PEllBig);
    PEllSplit := InvertibleSubmatrix(PEllBigSplit : epsinv := epsinv);
    // We create an approximately integral matrix from this; instead of
    // inverting we use NumericalLeftSolve.
    LMatAn := Matrix([ NumericalLeftSolve(PEllSplit, Matrix(row)) : row in Rows(PEllBigSplit) ]);
    LMat := Matrix([ [ FractionalApproximation(c : epscomp := epscomp, epsLLL := epsLLL) : c in Eltseq(row) ] : row in Rows(LMatAn) ]);
    M := Matrix(Basis(Lattice(LMat)));
    // Change to basis of lattice and return corresponding periods:
    Append(~Ls, Eltseq(CombinePeriodMatrix(M * PEllSplit)));
    Append(~col_numbers, col_number);
end for;

return Ls, col_numbers;

end function;
