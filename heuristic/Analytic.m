/***
 *  Finding an approximate basis for the geometric endomorphism ring
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


intrinsic ComplexStructure(P::.) -> .
{Gives the complex structure that corresponds to the period matrix P.}

CC := BaseRing(P); RR := RealField(CC);
PSplit := SplitPeriodMatrix(P); iPSplit := SplitPeriodMatrix(CC.1 * Matrix(k, P));
return NumericalLeftSolve(PSplit, iPSplit);

end intrinsic;


intrinsic RationalIsogenyEquations(JP::., JQ::.) -> .
{Given two complex structures JP and JQ, determines the corresponding equations satisfied by an isogeny between the two corresponding abelian varieties.}

// Basic invariants
RR := BaseRing(JP);
gP := #Rows(JP) div 2; gQ := #Rows(JQ) div 2;
n := 4 * gP * gQ;
// Building a formal matrix corresponding to all possible integral
// transformations of the lattice
R := PolynomialRing(CC, n);
vars := GeneratorsSequence(R);
Ma := Matrix(R, 2 * gP, 2 * gQ, vars);
// Condition that integral transformation preserve the complex structure
Comm := Eltseq(JP * Ma - Ma * JQ);
// Splitting previous linear equations by formal variable
return Matrix(RR, [ [MonomialCoefficient(c, var) : c in Comm] : var in vars ]);

end intrinsic;


intrinsic AnalyticRepresentationIsogeny(R::., P::., Q::.) -> .
{Given a rational representation R and two period matrices P and Q, finds an analytic representation of that same isogeny.}

// FIXME: We may not want to recalculate this every time and pass on P0 and s0
// as data. On the other hand, this is not a huge deal.
P0, s0 := InvertibleSubmatrix(P);
R := Matrix(BaseRing(P), R);
RowsRQ := Rows(R * Q);
RQ0 := Matrix(BaseRing(P), [ Eltseq(RowsRQ[i]) : i in s0 ]);
// Invert and return; transposes intervene because of right action
return Transpose(NumericalLeftSolve(Transpose(P0), Transpose(RQ0)));

end intrinsic;


intrinsic GeometricIsogenyApproximations(P::., Q::.) -> .
{Starting from period matrices P and Q, determines isogenies between the corresponding abelian varieties.}

// Setting up the equations
gP := #Rows(Transpose(P)); gQ := #Rows(Transpose(Q));
JP := ComplexStructure(P); JQ := ComplexStructure(Q);
// FIXME: Understand why this fails
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

As := [ AnalyticRepresentationIsogeny(R, P, Q) : R in Rs ];
return As, Rs;

end intrinsic;


intrinsic GeometricEndomorphismBasisApproximations(P::.) -> .
{Starting from a period matrix P, determines the endomorphisms of the corresponding abelian variety.}

return GeometricIsogenyApproximations(P, P);

end intrinsic;
