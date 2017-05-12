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
CC = ComplexField(prec*log(10)/log(2))
R.<x> = PolynomialRing(CC)

# Curve input: specify g and h in its equation y^2 + h y = g.
u = CC(-6).sqrt()

f = 1/7*(4*u - 60)*x^6 + 1/49*(-1536*u - 3840)*x^5 + 1/343*(-99840*u - 136704)*x^4 + 1/2401*(-1310720*u - 65536)*x^3 + 1/16807*(286457856*u + 85131264)*x^2 + 1/117649*(-1283457024*u + 2956984320)*x + 1/823543*(-2826960896*u - 3162505216)
f = 1/673*(-112*u - 6492)*x^6 + 1/452929*(5107200*u - 773952)*x^5 + 1/304821217*(19139022336*u + 2267408256)*x^4 + 1/205144679041*(878383923200*u + 17322438815744)*x^3 + 1/138062368994593*(-24289637111513088*u - 391596763165642752)*x^2 + 1/92915974333361089*(-21947333131409817600*u + 10080807600431235072)*x + 1/62532450726352012897*(-3864765776392296595456*u + 4319690447597534314496)

h = 0

# The main functionality:
End = EndomorphismData(f, h, prec = prec)
AsAlg, As, Rs = End.geometric_representations()
M = AsAlg[len(AsAlg)]
print End.field_of_definition()
print End.geometric().description()

#exit()

