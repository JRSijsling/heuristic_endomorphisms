/***
 *  Some LLL subtleties
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

// This load statement has to be changed:
//load "Algorithms.m";

// Testing third subalgorithm on a simple matrix
RR := RealField();
C := [
3.0001,4.9999,5.0001,
2.9999,2.0001,3.9999
];
M := Matrix(RR,2,3,C);

"Integral kernel approximation for a simple matrix:";
M;
for i:=0 to 8 do
    e := 2^i;
    "Exponent = ", e;
    L,K := IntegralKernelApproximation(M : N := 2^e);
    K;
end for;
