/***
 *  Polarizations and Rosati involutions
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function StandardSymplecticMatrix(g);

A := ScalarMatrix(g, 0);
B := ScalarMatrix(g, -1);
C := -B;
D := A;
return VerticalJoin(HorizontalJoin(A, B), HorizontalJoin(C, D));

end function;


function RosatiInvolution(AsAlg, As, Rs, A);
// Rosati involution

if IsExact(Parent(A)) then
    B := AsAlg;
    s := Eltseq(MatrixRatInBasisOverNF(A, B));
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

end function;


function DegreeEstimate(AsAlg, As, Rs, A);
// Round because of non-exact case
// Alternatively, could use Rdagger instead of Adagger

Adagger := RosatiInvolution(AsAlg, As, Rs, A);
tr := Trace(A * Adagger) * Factorial(#Rows(A) - 1);
if IsExact(Parent(A)) then
    return (Integers() ! tr);
else
    return Round(tr);
end if;

end function;
