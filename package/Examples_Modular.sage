"""
 *  Some examples of endomorphism rings of genus 2 hyperelliptic curves
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
"""

# Add if no initilization script set:
load('Initialize.sage')

# Ambient ring:
R.<x> = PolynomialRing(QQ)

# Curve input: specify g and h in its equation y^2 + h y = g.
# Hyperelliptic:
f = x^8 - 12*x^7 + 50*x^6 - 108*x^5 + 131*x^4 - 76*x^3 - 10*x^2 + 44*x - 19
h = 0

End = EndomorphismData(f, h, prec = 200)
AsAlg, As, Rs = End.geometric_representations()
print AsAlg
print Rs

Ds, idemsAlg, idemsAn = magma.SplittingInfoOneOff(AsAlg, As, Rs, nvals = 3)
Lats, col_number = magma.LatticesFromIdempotents([ idemsAn[1] ], End.period_matrix(), nvals = 2)
Lat = Lats[1]
print idemsAlg[1]
print col_number

gp.set_precision(prec)
g4 = gp.elleisnum(Lat.sage(), 4, flag = 1)
g6 = gp.elleisnum(Lat.sage(), 6, flag = 1)
eisvals_an = [ g4 , g6 ]
eisvals_alg = [ -26624/3, 741376/27 ]
CC = ComplexField()
print eisvals_an
print [ CC(eisval_alg) for eisval_alg in eisvals_alg ]

E = magma.EllipticCurve([ -eisvals_alg[0]/4, -eisvals_alg[1]/4 ])
print E
print Lat
print magma.Periods(E)

