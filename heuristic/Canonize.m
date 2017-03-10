/***
 *  Canonizing a matrix
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

function CanonizeMatrices(grep, frep, idems);

R<x> := PolynomialRing(Rationals());
f := R ! frep;
// TODO: Magma greatly loses here since the rationals do not admit a subfield constructor
if frep eq [0, 1] then
    L := Rationals();
    K := Rationals();
else
    L := NumberField(f);
    K := sub< L | L ! grep >;
end if;

return [ Matrix(K, [ [ K ! Eltseq(c) : c in Eltseq(row) ] : row in Rows(idem) ]) : idem in idems ];

end function;
