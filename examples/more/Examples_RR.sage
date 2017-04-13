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
a = -5/7
a = 4/3
# TODO: Bug
#a = 8/27
f = x^6 + 6*a^2*x^5 + 6*a^4*x^4 + (2*a^12 - 12*a^6 + 2)*x^3 + (3*a^14 - 6*a^8 + 3*a^2)*x^2 + (-6*a^16 + 12*a^10 - 6*a^4)*x + 2*a^18 - 4*a^12 + 2*a^6
f = 11*f
h = 0

CC.<I> = ComplexField(prec*log(10)/log(2))
R.<x> = PolynomialRing(CC)
r = CC(5).nth_root(3)
f = -16/27*x^6 - 32/9*r*x^5 - 32/9*r^2*x^4 + 1136/27*x^3 + 248/9*r*x^2 - 496/9*r^2*x + 2480/27
h = 0

# The main functionality:
End = EndomorphismData(f, h, prec = prec)
AsAlg, As, Rs = End.geometric_representations()
M = AsAlg[len(AsAlg)]

print End.geometric().description()

