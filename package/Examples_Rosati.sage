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

# LMFDB example 2:
f = 6*x^5 + 9*x^4 - x^3 - 3*x^2
h = 1

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
#K.<r> = NumberField(x^3 - x - 1)
#K.<r> = NumberField(x^2 - 5*x + 3)
#print End.over_field(K)
#print End.over_field(K).representations()
#print End.over_field(K).description()
print End.lattice()

# Verification of geometric endomorphisms:
print End.geometric_representations_check()

#exit()

