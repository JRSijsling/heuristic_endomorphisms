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
w = CC((1 + sqrt(-3))/2)
mu = 2 + w

f = x^6 + 6*mu^2*x^5 + 6*mu^4*x^4 + (2*mu^12 - 12*mu^6 + 2)*x^3 + (3*mu^14 -6*mu^8 + 3*mu^2)*x^2 + (-6*mu^16 + 12*mu^10 - 6*mu^4)*x + 2*mu^18 -4*mu^12 + 2*mu^6
h = 0

# The main functionality:
End = EndomorphismData(f, h, prec = prec)
AsAlg, As, Rs = End.geometric_representations()
M = AsAlg[len(AsAlg)]
print End.field_of_definition()
print End.geometric().description()

#exit()

