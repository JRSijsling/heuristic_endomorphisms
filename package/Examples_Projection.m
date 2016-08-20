/***
 *  Examples of the projection functionality
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

load "Projection.m";

R<x> := PolynomialRing(Rationals());

fX := x^6 + x^4 + x^2 + 1;
fE := x^3 + x^2 + x + 1;
A := [0, 2];
deg := 2;

fX := x^6 + x^4 + x^2 + 1;
fE := x^4 + x^3 + x^2 + x;
A := [0, 2];
deg := 2;

fX := x^5 + x^3 + x;
fE := x*(x - 2)*(x + 1);
A := [-1, 1];
deg := 2;

fX := x^5 + x^3 + x;
fE := (x + 1)*(x - 1)*(x + 2);
A := [-1, 1];
deg := 2;

fX := x^5 + x^3 + x;
fE := -2*x^4 - x^3 + 2*x^2 + x;
A := [1, -1];
deg := 2;

fX := x^6 - 8*x^4 + 2*x^3 + 16*x^2 - 36*x - 55;
fE := -10582/27*x^4 + 215/3*x^3 + x;
A := [1, 1];
deg := 2;

//time MyTest();
time test, proj := ProjectionToEllipticFactorG2(fX, fE, A, deg);
print proj;

exit;
