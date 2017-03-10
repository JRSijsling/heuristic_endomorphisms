/***
 *  PEL structures: this functionality still has to be added and is now out of
 *  date with the main algorithms (such as the use of IntegralKernel)
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

function TransformPELData(E,T);
// As of yet, this only involves the polarization

if IsZero(E) then
    F := 0;
else
    // We have an integral matrix of unit determinant, so considering that Magma
    // preserves base ring, we can invert without fear of instability.
    Ti := T^(-1);
    F := Transpose(Ti)*E*Ti;
end if;

return F;

end function;


function NeronSeveriEquations(J);

// Basic invariants
k := Parent(J[1,1]);
g := Integers()!(#Rows(J)/2);
n := 4*g^2;

// Building a formal matrix corresponding to all possible integral
// transformations of the lattice
R := PolynomialRing(k,n);
Ca := [];
for i:=1 to n do
    Append(~Ca,R.i);
end for;
Ma := Matrix(R,2*g,2*g,Ca);

// Alternating condition
D1 := Ma + Transpose(Ma);
// Hermitian condition
D2 := Transpose(J)*Ma + Ma*J;

//Resulting linear equations
C := [];
for d in Eltseq(D1) cat Eltseq(D2) do
    for i:=1 to n do
        Append(~C,MonomialCoefficient(d,R.i));
    end for;
end for;

// Constructing and returning matrix
M := Matrix(k,2*n,n,C);
return M;

end function;


function NeronSeveriBasis(P0 : eps := 2^(-32));

// Basic invariants
k := Parent(P0[1,1]);
g := #Rows(Transpose(P0));

// Optimize
P,T := OptimizedPeriodMatrix(P0);
P := P0; T := IdentityMatrix(Integers(),2*g);

// Setting up the equations
J := ComplexStructure(P);
M := NeronSeveriEquations(J);
L,K := IntegralKernel(M);

Es := [];
for r in Rows(K) do
    E := Matrix(Integers(),2*g,2*g,Eltseq(r));
    // Alternating condition
    D1 := E + Transpose(E);
    // Hermitian condition
    D2 := Transpose(J)*E + E*J;
    C := Eltseq(D1) cat Eltseq(D2);
    // Culling the correct transformations
    if &and([Abs(c) lt eps : c in C]) then
        Append(~Es,E);
    end if;
end for;

// Transform back from optimization
Es := [Transpose(T)*E*T : E in Es];
return Es;

end function;


function Polarizations(P);
// Add positivity condition

return [0];

end function;


