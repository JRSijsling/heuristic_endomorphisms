/***
 *  Recognizing complex numbers as algebraic numbers in relative fields
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


intrinsic RelativeMinimalPolynomial(a::FldComElt, F::FldNum) -> RngUPolElt
{Determines a relative minimal polynomial of the element a with respect to the stored infinite place of F.}

// NOTE: The speed of this relies on Magma storing the evaluation of a
// generator (called iota below). It seems to do this well. Otherwise we should
// store it as part of the structure of F.
degF := Degree(F); R<x> := PolynomialRing(F);
CC := Parent(a); RR := RealField(CC); prec := Precision(CC);

// NOTE: Here height is a parameter to play with.
degf := 0; height := 1; height0 := 100;
gen := CC ! Evaluate(F`iota, F.1 : Precision := prec);
powersgen := [ gen^i : i in [0..(degF - 1)] ];
powersa := [ CC ! 1 ];
MLine := [ ];

while true do
    // Increase height and number of possible relations
    degf +:= 1; height *:= height0;
    powera := powersa[#powersa] * a; powersa cat:= [ powera ];
    MLine cat:= [ powersgen * a : powergen in powersgen ];
    M := Transpose(Matrix(CC, [ MLine ]));

    // Now split and take an IntegralLeftKernel
    MSplit := SplitMatrix(M);
    Ker := IntegralLeftKernel(MSplit);
    for row in Rows(Ker) do
        // Height condition
        if &and[ Abs(c) le height : c in Eltseq(row) ] then
            f := &+[ &+[ row[i*degF + j + 1]*F.1^j : j in [0..(degF - 1)] ] * x^i : i in [0..degf] ];
            // Factor (to eliminate redundancy) and check
            Fac := Factorization(f);
            for tup in Fac do
                g := tup[1];
                if (RR ! Abs(Evaluate(g, a))) lt RR`epscomp then
                    return g;
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

CC := Parent(a); RR := RealField(CC);
M := Matrix(RR, [ [ 1 ], [ -Real(a) ] ]);
K := IntegralLeftKernel(M); q := K[1,1] / K[1,2];
if (RR ! Abs(q - a)) lt RR`epscomp then
    return q;
else
    error Error("LLL not does return a sufficiently good fraction");
end if;

end intrinsic;


intrinsic AlgebraizeElementInRelativeField(a::FldComElt, K::FldNum) -> FldNumElt
{Finds an algebraic approximation of a as an element of K.}

degK := Degree(K); R<x> := PolynomialRing(K);
F := BaseField(K); degF := Degree(K);
CC := Parent(a); RR := RealField(CC); prec := Precision(CC);

genK := CC ! Evaluate(K`iota, K.1 : Precision := prec); genF := CC ! Evaluate(F`iota, F.1 : Precision := prec);
powersgenK := [ genK^i : i in [0..(degK - 1)] ]; powersgenF := [ genF^i : i in [0..(degF - 1)] ];
MLine := &cat[ [ powergenF * powergenK : powergenF in powersgenF ] : powergenK in powersgenK ] cat [-a];
M := Transpose(Matrix(CC, [ MLine ]));

// Now split and take an IntegralLeftKernel
MSplit := SplitMatrix(M);
Ker := IntegralLeftKernel(MSplit);
for row in Rows(Ker) do
    den := row[#row];
    if den ne 0 then
        sCC := &+[ &+[ row[i*degF + j]*genF^j : j in [0..(degF - 1)] ] * genK^i : i in [0..(degK - 1)] ] / den;
        // Check correct to given precision
        if (RR ! Abs(sCC - a)) lt RR`epscomp then
            s := &+[ &+[ row[i*degF + j]*F.1^j : j in [0..(degF - 1)] ] * K.1^i : i in [0..(degK - 1)] ] / den;
            return s;
        end if;
    end if;
end for;
error Error("LLL fails to algebraize element in ambient");

end intrinsic;


intrinsic AlgebraizeMatrixInRelativeField(A::AlgMatElt, K::FldNum) -> AlgMatElt
{Algebraizes a matrix.}

return Matrix([ [ AlgebraizeElementInRelativeField(c, K) : c in Eltseq(row) ] : row in Rows(A) ]);

end intrinsic;


intrinsic AlgebraizeMatrixInRelativeField(A::ModMatRngElt, K::FldNum) -> ModMatRngElt
{Algebraizes a matrix.}

return Matrix([ [ AlgebraizeElementInRelativeField(c, K) : c in Eltseq(row) ] : row in Rows(A) ]);

end intrinsic;


intrinsic GeometricEndomorphismBasisRepresentations(As::SeqEnum, Rs::SeqEnum, K::FldNum) -> List
{Final algebraization step.}

AsAlg := [ AlgebraizeMatrixInRelativeField(A, K) : A in As ];
return [* As, Rs, AsAlg *];

end intrinsic;
