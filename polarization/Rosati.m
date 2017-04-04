/***
 *  Polarizations and Rosati involutions
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


intrinsic StandardSymplecticMatrix(g::RngIntElt) -> .
{Standard symplectic 2 g x 2 g matrix.}

A := ScalarMatrix(g, 0); B := ScalarMatrix(g, -1); C := -B; D := A;
return VerticalJoin(HorizontalJoin(A, B), HorizontalJoin(C, D));

end intrinsic;


intrinsic RosatiInvolution(GeoEndList::List, A::.) -> .
{Returns the Rosati involution of A.}

AsAlg, Rs, As := Explode(GeoEndList);
if IsExact(Parent(A)) then
    B := AsAlg;
    s := Eltseq(MatrixInBasis(A, B));
else
    B := As;
    s := Eltseq(MatrixInBasis(A, B));
    s := [ Round(c) : c in s ];
end if;
R := &+[ s[i] * Rs[i] : i in [1..#Rs] ];
J := StandardSymplecticMatrix(#Rows(A));
Rdagger := -J * Transpose(R) * J;
sdagger := Eltseq(MatrixInBasis(Rdagger, Rs));
Adagger := &+[ sdagger[i] * B[i] : i in [1..#Rs] ];
return Adagger;

end intrinsic;


intrinsic DegreeEstimate(GeoEndList::List, A::.) -> .
{Estimates degree of corresponding endomorphism.}

Adagger := RosatiInvolution(GeoEndList, A);
tr := Trace(A * Adagger) * Factorial(#Rows(A) - 1);
if IsExact(Parent(A)) then
    return (Integers() ! tr);
else
    return Round(tr);
end if;

end intrinsic;
