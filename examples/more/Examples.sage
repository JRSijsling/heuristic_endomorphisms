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

# Some interesting curves:
# Endomorphisms and decompositions defined over QQ:
f = x^4 + x^2
h = x^3 + 1
# Case where [RR, CC] occurs:
f = x^5 + x^4 + 2*x^3 + x^2 + x
h = x^2 + x
# The big degree 48 case:
f = x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1
h = R(0)
# A case with trivial endomorphism ring that shows an alternative input method:
f = [-8,12,8,-8,-8,1,1]
h = [0]
# Non-cyclic CM:
f = x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8
h = 0
# RM over QQ:
f = -x^5
h = x^3 + x + 1
# Potential RM:
f = x^6 + 2*x^3 - x
h = x^3 + 1
# Debugging the splitting functionality:
f = 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x
h = x
# LMFDB example 1:
f = x^6 + 2*x^3 - x
h = x^3 + 1
# LMFDB example 2:
f = 6*x^5 + 9*x^4 - x^3 - 3*x^2
h = 1
# LMFDB example 3:
f = -2*x^4 + 4*x^2 - 9*x - 14
h = x^3 + 1
# Debug (TODO: FindPoint):
f = 3*x^3 - 2*x^2 + 6*x + 2
h = x^3 + x
# Debug (TODO: Subfield):
f = x^4 + x^3 + 2*x^2 + x + 1
h = x^3 + x^2 + x + 1
# Debug:
f = x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8
h = 0

# Apply a substition if necessary to get around Magma bugs (commented out for now)
#f = R(f)
#h = R(h)
#subst = x + 1
#subst = -2/3*x + 1
#subst = (-1/3*x + 1/3)/(-3/2*x + 2)
#den = subst.denominator()
#f = R(den^6 * f(subst))
#h = R(den^3 * h(subst))

# The main functionality:
End = EndomorphismData(f, h, prec = 200)
AsAlg, As, Rs = End.geometric_representations()
M = AsAlg[len(AsAlg)]
#print End
#print M
#print End.rosati_involution(M)
#print End.degree_estimate(M)
print End.field_of_definition()
#print End.geometric()
print End.geometric().representations()[0]
print End.geometric().description()
#print End.over_base()
#print End.over_base().representations()
#print End.over_base().description()
#K.<r> = NumberField(x^2 - 5*x + 3)
#K.<r> = NumberField(x^2 - x + 1)
#print End.over_field(K)
#print End.over_field(K).representations()
#print End.over_field(K).description()
print End.lattice()

# Verification of decomposition data:
#Dec = End.decomposition()
#print Dec
#print Dec.field_of_definition()
#print Dec._idems_[1]
#print Dec.factors()
#print Dec.certificate_g2()

# Verification of geometric endomorphisms:
print End.geometric_representations_check()

#exit()

