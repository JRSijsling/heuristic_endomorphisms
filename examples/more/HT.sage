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
t = CC(5/2)
s = CC(-3 + 14*t^2 - 27*t^4).sqrt()

A = s / (2*t)
B = (1 + 3*t^2)/(1 - 3*t^2)
Q = -(1 - 2*t^2 + 9*t^4)*(1 - 28*t^2 + 166*t^4 - 252*t^6 + 81*t^8)/(4*t^2*(1 - 3*t^2)^2*(1 - t^2)*(1 - 9*t^2))
f = x*(x^4 + (A - B)*x^3 + Q*x^2 + (A + B)*x + 1)
h = 0

# The main functionality:
End = EndomorphismData(f, h, prec = prec)
AsAlg, As, Rs = End.geometric_representations()
M = AsAlg[len(AsAlg)]
print End.field_of_definition()
print End.geometric().description()

#exit()

