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

load('Initialize.sage')

# Ambient ring:
R.<x> = PolynomialRing(QQ)

# Curve input: specify g and h in its equation y^2 + h y = g.

# Some interesting curves:
# Many automorphisms, and descent to a smaller field down from degree 12:
# (a personal favorite)
g = R(1)
h = x^3 + 1
# Another degree 12 field:
g = -2*x^6 + 2*x^3 - 1
h = x^3 + 1
# Endomorphisms and decompositions defined over QQ:
g = x^4 + x^2
h = x^3 + 1
# Case where [RR, CC] occurs:
g = x^5 + x^4 + 2*x^3 + x^2 + x
h = x^2 + x
# The big degree 48 case:
#g = x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1
#h = R(0)
# Non-cyclic CM:
g = x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8
h = 0
# A case with trivial endomorphism ring that shows an alternative input method:
#g = [-8,12,8,-8,-8,1,1]
#h = [0]

# Apply a substition if necessary to get around Magma bugs (commented out for now)
#g = R(g)
#h = R(h)
#subst = x + 1
#subst = -2/3*x + 1 
#subst = (-1/3*x + 1/3)/(-3/2*x + 2)
#den = subst.denominator()
#g = R(den^6 * g(subst))
#h = R(den^3 * h(subst))

# The main functionality:
End = EndomorphismData(g, h, prec = 200)
K.<r> = NumberField(x^3 - x - 1)
K.<r> = NumberField(x^2 - 5*x + 3)
print End
print End.geometric_representations()
print End.field_of_definition()
print End.geometric()
print End.geometric().representations()
print End.geometric().description()
print End.over_base()
print End.over_base().representations()
print End.over_base().description()
print End.over_field(K)
print End.over_field(K).representations()
print End.over_field(K).description()
print End.lattice()
print End.decomposition()
print End.decomposition().field_of_definition()
print End.decomposition().factors()
