/***
 *  Finding an approximate basis for the geometric endomorphism ring through LLL
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

function ComplexStructure(P);
// Input:   A period matrix of dimension 2g x g.
// Output:  A 2g by 2g matrix left multiplication by which represents
//          multiplication by i on the basis furnished by the columns.

// Basic invariants
k := BaseRing(P);
i := k.1;

// Splitting both P and i*P
PSplit := SplitPeriodMatrix(P);
iPSplit := SplitPeriodMatrix(i * Matrix(k, P));
J := NumericalSolve(PSplit, iPSplit);

return J;

end function;


function RationalEndomorphismEquations(J);
// Input:   A matrix J over RR of dimension 2g x 2g, representing multiplication
//          by i.
// Output:  The corresponding equations for the entries of a commuting matrix in
//          the rational representation, as a matrix of which we take the left
//          kernel; later we take integral such solutions, that like J will act
//          on the left.

// Basic invariants
k := BaseRing(J);
g := #Rows(J) div 2;
n := 4 * g^2;

// Building a formal matrix corresponding to all possible integral
// transformations of the lattice
R := PolynomialRing(k, n);
vars := GeneratorsSequence(R);
Ma := Matrix(R, 2 * g, 2 * g, vars);

// Condition that integral transformation preserve the complex structure
Comm := Eltseq(J * Ma - Ma * J);

// Integration of PEL algorithms; forget about this for now
//if not IsZero(E) then
//    Comm := Comm cat Eltseq(Transpose(Ma) * E - E * Ma);
//end if;

// Splitting previous linear equations by formal variable
M := Matrix([[MonomialCoefficient(c, var) : c in Comm] : var in vars]);

// Order: in the end, we want solutions (a, b, c, d) (or larger) in a left
// kernel. For this to make sense the first rows has to correspond with a, et
// cetera.

return M;

end function;


function AnalyticRepresentation(R, P : J := 0, P0, s0, epsinv := epsinv0);
// Input:   A matrix R giving the rational representation of an endomorphism,
//          along with the period matrix P. Extra arguments speed up the
//          computation.
// Output:  The analytic representation A of the corresponding endomorphism,
//          acting on the right. So R * P = P * A.

// Calculate additional data if necessary
if IsZero(J) then
    J := ComplexStructure(P);
    P0, s0 := InvertibleSubmatrix(P : epsinv := epsinv);
end if;

// Passing from rational to analytic by the invertible submatrix;
// this could be improved by taking the most stable of these
R := Matrix(BaseRing(P), R);
RowsRP := Rows(R * P);
RP0 := Matrix([RowsRP[i] : i in s0]);

// Invert and return; transposes intervene because of right action
return Transpose(NumericalSolve(Transpose(P0), Transpose(RP0)));

end function;


function PeriodMatrix(h, g : prec := prec0)
// Input:   Two polynomials h and g over the rationals, reals or CC.
// Output:  The corresponding period matrix.

R := PolynomialRing(ComplexField(prec));
f := (R!g) + ((R!h)/2)^2;
J := AnalyticJacobian(f);
return Transpose(BigPeriodMatrix(J));

end function;


function GeometricEndomorphismBasisFromPeriodMatrix(P : epscomp := epscomp0,
    epsLLL := epsLLL0, epsinv := epsinv0)
// Input:   A period matrix of dimension 2g x g.
// Output:  An approximate basis of the corresponding endomorphism ring,
//          returned in both analytic and rational representations.

// Optimization and PEL data used to be inserted here; add later
//P, T := OptimizedPeriodMatrix(P0);
//E := TransformPELData(E0, T);

// Setting up the equations
g := #Rows(Transpose(P));
J := ComplexStructure(P);
P0, s0 := InvertibleSubmatrix(P : epsinv := epsinv);
M := RationalEndomorphismEquations(J);

// Determination of approximate endomorphisms by LLL
K := IntegralKernel(M : epsLLL := epsLLL);

// Deciding which rows to keep
Rs := [];
for r in Rows(K) do
    alpha := Matrix(Integers(), 2*g, 2*g, Eltseq(r));

    // Another PEL point
    //if not IsZero(E) then
    //    D2 := Transpose(alpha)*E - E*alpha;
    //    C := C cat Eltseq(D2);
    //end if;

    // Culling the correct transformations from holomorphy condition
    Comm := alpha * J - J * alpha;
    if &and([Abs(c) lt epscomp : c in Eltseq(Comm)]) then
        Append(~Rs, alpha);
    end if;
end for;

// PEL: Transform back from optimization
//Rs := [T^(-1)*R*T : R in Rs];
// The following line is a sanity check and should cause an incorrect response
// for the first entry in the LMFDB, causing the algebraization of the matrices
// to take too long.
//Rs := [Transpose(R) : R in Rs];
As := [AnalyticRepresentation(R, P : J := J, P0 := P0, s0 := s0, 
    epsinv := epsinv) : R in Rs];

return As, Rs;

end function;
