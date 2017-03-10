/***
 *  Polarizations and Rosati involutions
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
