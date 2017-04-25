/***
 *  Finding an approximate basis for the geometric endomorphism ring
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic ComplexStructure(P::.) -> .
{Gives the complex structure that corresponds to the period matrix P.}

CC := BaseRing(P); RR := RealField(CC);
PSplit := SplitMatrix(P); iPSplit := SplitMatrix(CC.1 * P);
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
R := PolynomialRing(RR, n);
vars := GeneratorsSequence(R);
M := Matrix(R, 2 * gP, 2 * gQ, vars);
// Condition that integral transformation preserve the complex structure
Comm := Eltseq(ChangeRing(JP, R) * M - M * ChangeRing(JQ, R));
// Splitting previous linear equations by formal variable
return Matrix(RR, [ [MonomialCoefficient(c, var) : c in Comm] : var in vars ]);

end intrinsic;


intrinsic AnalyticRepresentationIsogeny(R::., P::., Q::.) -> .
{Given a rational representation R of an isogeny and two period matrices P and Q, finds an analytic representation of that same isogeny.}

// FIXME: We may not want to recalculate this every time and pass on P0 and s0
// as data. On the other hand, this is not a huge deal.
P0, s0 := InvertibleSubmatrix(P);

R := Matrix(BaseRing(P), R);
RowsRQ := Rows(R * Q);
RQ0 := Matrix(BaseRing(P), [ Eltseq(RowsRQ[i]) : i in s0 ]);
// Invert and return; transposes intervene because of right action
return Transpose(NumericalLeftSolve(Transpose(P0), Transpose(RQ0)));

end intrinsic;


intrinsic AnalyticRepresentationEndomorphism(R::., P::.) -> .
{Given a rational representation R of an endomorphism and a period matrix P, finds an analytic representation of that same endomorphism.}

return AnalyticRepresentationIsogeny(R, P, P);

end intrinsic;


intrinsic GeometricIsogenyBasisApproximations(P::., Q::.) -> .
{Starting from period matrices P and Q, determines isogenies between the corresponding abelian varieties.}

// Basic invariants
gP := #Rows(Transpose(P)); gQ := #Rows(Transpose(Q));
JP := ComplexStructure(P); JQ := ComplexStructure(Q);

// FIXME: Understand why this fails
//JP := Transpose(JP); JQ := Transpose(JQ);
M := RationalIsogenyEquations(JP, JQ);

// Determination of approximate endomorphisms by LLL
Ker := IntegralLeftKernel(M);

// Deciding which rows to keep
RR := BaseRing(JP); gensHom := [];
for r in Rows(Ker) do
    genHom := Matrix(Rationals(), 2*gP, 2*gQ, Eltseq(r));
    // Culling the correct transformations from holomorphy condition
    Comm := JP * ChangeRing(genHom, RR) - ChangeRing(genHom, RR) * JQ;
    if &and([ (RR ! Abs(c)) lt RR`epscomp : c in Eltseq(Comm) ]) then
        Append(~gensHom, genHom);
    end if;
end for;

gensTan := [ AnalyticRepresentationEndomorphism(genHom, P) : genHom in gensHom ];
return [* gensTan, gensHom *];

end intrinsic;


intrinsic GeometricEndomorphismBasisApproximations(P::.) -> .
{Starting from a period matrix P, determines the endomorphisms of the corresponding abelian variety.}

return GeometricIsogenyBasisApproximations(P, P);

end intrinsic;
