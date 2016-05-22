/***
 *  Elliptic curves from decompositions
 *
 *  Copyright (C) 2016  J.R. Sijsling (sijsling@gmail.com)
 *
 *  Distributed under the terms of the GNU General License (GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc., 51
 *  Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

function SmallIdempotent(Q, GensQ);
// Input:   A split quaterion algebra and a basis for an order.
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

Ls := [];
for idem in idemsAn do
    // Create analytic idempotent and project:
    PEllHuge := P * idem;
    PEllBig := Transpose(Matrix([ [ row[1] : row in Rows(PEllHuge) ] ]));
    // Change row if the result is too small (this happens!):
    if &and([ Abs(c) lt epscomp : c in Eltseq(PEllBig) ]) then
        PEllBig := Transpose(Matrix([ [ row[2] : row in Rows(PEllHuge) ] ]));
    end if;
    PEllBigSplit := SplitPeriodMatrix(PEllBig);
    PEllSplit := InvertibleSubmatrix(PEllBigSplit : epsinv := epsinv);
    // We create an approximately integral matrix from this; instead of
    // inverting we use NumericalSolve.
    LMatAn := Matrix([ NumericalSolve(PEllSplit, Matrix(row)) : row in
        Rows(PEllBigSplit) ]);
    LMat := Matrix([ [ FractionalApproximation(c : epscomp := epscomp, epsLLL :=
        epsLLL) : c in Eltseq(row) ] : row in Rows(LMatAn) ]);
    M := Matrix(Basis(Lattice(LMat)));
    // Change to basis of lattice and return corresponding periods:
    Append(~Ls, Eltseq(CombinePeriodMatrix(M * PEllSplit)));
end for;

return Ls;

end function;


// Deprecated:
function EllipticCurveFromEisensteinValues(eisvalsalg);
// Input:   An algebraized tuple of outputs of elleisnum.
// Output:  The corresponding elliptic curve.

if #eisvalsalg[1] eq 1 then
    // Rational case: Cremona label.
    E := EllipticCurve([-eisvalsalg[1][1]/48, -eisvalsalg[2][1]/864]);
    return CremonaReference(E);
else
    // Otherwise return the Eltseq as a Weierstrass equation:
    // to be divided by -48 and -864.
    return eisvalsalg;
end if;

end function;

// TODO: Recover certain invariants like field of definition without having to
// go through the whole EndomorphismLattice business, and make this a function
// in its own right.
