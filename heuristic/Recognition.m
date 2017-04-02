/***
 *  Functionality for recognizing complex numbers as algebraic numbers and
 *  related optimizations
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


intrinsic RelativeMinimalPolynomial(a::FldComElt, F::FldNum, iota::PlcNumElt) -> RngUPolElt
{Determines a relative minimal polynomial of the element a with respected to the number field F and the embedding iota of F.}
// TODO: Keep track of loop by a verbose flag

CC := Parent(a); RR := RealField(CC); prec := Precision(CC);
epscomp := 10^(-prec + 30); epsLLL := 5^(-prec + 7);
R<x> := PolynomialRing(F);

// NOTE: h0 is a parameter to play with
d := 0; h := 1; h0 := 10;
powers := [ CC ! 1 ];
while true do
    // Increase height:
    d +:= 1; h *:= h0;
    Append(~powers, powers[#powers] * a);
    M := Transpose(Matrix([powers]));
    // Now split and take an IntegralLeftKernel:
    MSplit := SplitPeriodMatrix(M);
    Ker := IntegralLeftKernel(MSplit : epsLLL := epsLLL);
    for row in Rows(Ker) do
        // Height condition:
        if &and[ Abs(c) le h : c in Eltseq(row) ] then
            pa := &+ [ row[i + 1] * x^i : i in [0..d] ];
            // Factor (to eliminate redundancy) and check:
            Fac := Factorization(pa);
            for tup in Fac do
                fa := tup[1];
                if Abs(Evaluate(fa, a)) lt epscomp then
                    return fa;
                end if;
            end for;
        end if;
    end for;
end while;

end intrinsic;


intrinsic RelativeMinimalPolynomials(L::SeqEnum) -> SeqEnum
{Polynomializes matrices.}

return [ RelativeMinimalPolynomial(a) : a in L ] ;

end intrinsic;


intrinsic RelativeMinimalPolynomialsMatrices(As::SeqEnum) -> SeqEnum
{Polynomializes matrices.}

return &cat[ RelativeMinimalPolynomials(Eltseq(A)) : A in As ];

end intrinsic;


intrinsic FractionalApproximation(a::FldComElt) -> FldRatElt
{Fractional approximation of a complex number a.}
    CC := Parent(a); prec := Precision(CC);
    epscomp := 10^(-prec + 30); epsLLL := 5^(-prec + 7);
    M := Matrix([[1], [-Real(a)]]);
    K := IntegralLeftKernel(M : epsLLL := epsLLL);
    q := K[1,1] / K[1,2];
    if Abs(q - a) lt epscomp then
        return q;
    else
        error Error("LLL not does return a sufficiently good fraction");
    end if;
end intrinsic;







intrinsic AlgebraizeElementInRelativeField(a::FldComElt, K::FldNum, iota::PlcNumElt) -> FldNumElt
{Finds an algebraic approximation of a as an element of K via the embedding iota of K into CC.}

CC := Parent(a); RR := RealField(CC); prec := Precision(CC);
epscomp := 10^(-prec + 30); epsLLL := 5^(-prec + 7);
R<x> := PolynomialRing(K);

Ca := Parent(a);
ran := Ca!h;
Ra := RealField(Precision(Ca));
K := NumberField(Polynomial(frep));
d := Degree(K);
if d ne 1 then
    r := K.1;
else
    r := K ! 0;
end if;
// Creating algebraic and analytic powers:
powers := [ r^n : n in [0..d-1] ];
powersan := [ ran^n : n in [0..d-1] ];
// Column matrix of embedding of basis elements plus (minus) a:
Man := VerticalJoin(Transpose(Matrix([powersan])), Matrix([[-a]]));
// Now split and take an IntegralLeftKernel:
ManSplit := SplitPeriodMatrix(Man);
Ker := IntegralLeftKernel(ManSplit : epsLLL := epsLLL);
for R in Rows(Ker) do
    den := R[d + 1];
    if den ne 0 then
        s := (&+[ R[i] * powers[i] : i in [1..d] ])/ den;
        san := (&+[ R[i] * powersan[i] : i in [1..d] ]) / den;
        // Check correct to given precision:
        if Abs(san - a) lt epscomp then
            return Eltseq(s);
        end if;
    end if;
end for;

error Error("Failed to algebraize element in ambient");

end intrinsic;


intrinsic AlgebraizeMatrixInRelativeField(A::AlgMatElt, K::FldNum, iota::PlcNumElt) -> AlgMatElt
{Algebraizes a matrix.}

return Matrix([ [ AlgebraizeElementInRelativeField(c, K, iota) : c in Eltseq(row) ] : row in Rows(A) ]);

end intrinsic;


intrinsic AlgebraizeMatrixInRelativeField(A::ModMatRngElt, K::FldNum, iota::PlcNumElt) -> ModMatRngElt
{Algebraizes a matrix.}

return Matrix([ [ AlgebraizeElementInRelativeField(c, K, iota) : c in Eltseq(row) ] : row in Rows(A) ]);

end intrinsic;


intrinsic GeometricEndomorphismRepresentations(As::SeqEnum, Rs::SeqEnum, K::FldNum, iota::PlcNumElt) -> List
{Final algebraization step.}

AsAlg := [ AlgebraizeMatrixInRelativeField(A, K, iota) : A in As ];
return [* As, Rs, AsAlg *];

end intrinsic;
