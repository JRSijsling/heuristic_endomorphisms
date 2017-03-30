/***
 *  Finding an approximate basis for the geometric endomorphism ring through LLL
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
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
J := NumericalLeftSolve(PSplit, iPSplit);

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

// Splitting previous linear equations by formal variable
M := Matrix([[MonomialCoefficient(c, var) : c in Comm] : var in vars]);

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
return Transpose(NumericalLeftSolve(Transpose(P0), Transpose(RP0)));

end function;


function GeometricEndomorphismBasisFromPeriodMatrix(P : epscomp := epscomp0,
    epsLLL := epsLLL0, epsinv := epsinv0)
// Input:   A period matrix of dimension 2g x g.
// Output:  An approximate basis of the corresponding endomorphism ring,
//          returned in both analytic and rational representations.

// Setting up the equations
g := #Rows(Transpose(P));
J := ComplexStructure(P);
P0, s0 := InvertibleSubmatrix(P : epsinv := epsinv);
M := RationalEndomorphismEquations(J);

// Determination of approximate endomorphisms by LLL
K := IntegralLeftKernel(M : epsLLL := epsLLL);

// Deciding which rows to keep
Rs := [];
for r in Rows(K) do
    alpha := Matrix(Rationals(), 2*g, 2*g, Eltseq(r));

    // Culling the correct transformations from holomorphy condition
    Comm := alpha * J - J * alpha;
    if &and([Abs(c) lt epscomp : c in Eltseq(Comm)]) then
        Append(~Rs, alpha);
    end if;
end for;

As := [AnalyticRepresentation(R, P : J := J, P0 := P0, s0 := s0,
    epsinv := epsinv) : R in Rs];

// The end result is the actions As on the tangent space,
// which are the duals of the actions on the differentials.
// These are represented as a right multiplication because I see vectors as
// rows, just like Magma and Sage.
// The action the Rs is on the right, and in the end we always have
// R * P = P * A .
return As, Rs;

end function;


function RationalIsogenyEquations(JP, JQ);

// Basic invariants
k := BaseRing(JP);
gP := #Rows(JP) div 2; gQ := #Rows(JQ) div 2;
n := 4 * gP * gQ;

// Building a formal matrix corresponding to all possible integral
// transformations of the lattice
R := PolynomialRing(k, n);
vars := GeneratorsSequence(R);
Ma := Matrix(R, 2 * gP, 2 * gQ, vars);

// Condition that integral transformation preserve the complex structure
Comm := Eltseq(JP * Ma - Ma * JQ);

// Splitting previous linear equations by formal variable
M := Matrix([[MonomialCoefficient(c, var) : c in Comm] : var in vars]);

return M;

end function;


function AnalyticRepresentationIsogeny(R, P, Q : epsinv := epsinv0);

// Passing from rational to analytic by the invertible submatrix;
// this could be improved by taking the most stable of these
P0, s0 := InvertibleSubmatrix(P : epsinv := epsinv);
R := Matrix(BaseRing(P), R);
RowsRQ := Rows(R * Q);
RQ0 := Matrix(BaseRing(P), [ Eltseq(RowsRQ[i]) : i in s0 ]);

// Invert and return; transposes intervene because of right action
return Transpose(NumericalLeftSolve(Transpose(P0), Transpose(RQ0)));

end function;


function GeometricIsogenyBasisFromPeriodMatrices(P, Q : epscomp := epscomp0,
    epsLLL := epsLLL0, epsinv := epsinv0)

// Setting up the equations
gP := #Rows(Transpose(P)); gQ := #Rows(Transpose(Q));
JP := ComplexStructure(P); JQ := ComplexStructure(Q);
//JP := Transpose(JP); JQ := Transpose(JQ);
M := RationalIsogenyEquations(JP, JQ);

// Determination of approximate endomorphisms by LLL
K := IntegralLeftKernel(M : epsLLL := epsLLL);

// Deciding which rows to keep
Rs := [];
for r in Rows(K) do
    alpha := Matrix(Rationals(), 2*gQ, 2*gP, Eltseq(r));
    // Culling the correct transformations from holomorphy condition
    Comm := JP * alpha - alpha * JQ;
    if &and([Abs(c) lt epscomp : c in Eltseq(Comm)]) then
        Append(~Rs, alpha);
    end if;
end for;

As := [AnalyticRepresentationIsogeny(R, P, Q) : R in Rs];
return As, Rs;

end function;
