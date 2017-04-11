"""
 *  Some examples of endomorphism rings
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

# Ambient ring:
R.<x> = PolynomialRing(QQ)
#F.<r> = NumberField(x^2 - 2)
#R.<x> = PolynomialRing(F)

# Curve input: specify g and h in its equation y^2 + h y = g.

# The big degree 48 case:
f = x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1
h = R(0)
# Debug RM over QQ:
f = -x^5
h = x^3 + x + 1
# Debug splitting functionality:
f = 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x
h = x
# Debug (TODO: FindPoint):
#f = 3*x^3 - 2*x^2 + 6*x + 2
#h = x^3 + x
# Debug (TODO: Subfield):
#f = x^4 + x^3 + 2*x^2 + x + 1
#h = x^3 + x^2 + x + 1
# Debug (TODO: Point returned is wrong):
#f = x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8
#h = 0
# Genus 3:
f = x^7 + x^6 + x^5 + x^3 + x^2 + x
h = x^4 + x^2 + 1
# Case where [RR, CC] occurs:
f = x^5 + x^4 + 2*x^3 + x^2 + x
h = x^2 + x

X = HyperellipticCurve(f, h)

## Ambient ring:
#P2.<x,y,z> = ProjectiveSpace(QQ, 2)
#F = x^4 - 3*x^3*y + 3*x^3*z + 3*x^2*y^2 + 4*x^2*y*z + x^2*z^2 + 2*x*y^3 + 4*x*y^2*z + x*y*z^2 - x*z^3 - 2*y^4 + 4*y^3*z - 2*y^2*z^2 + 3*y*^3 - 2*z^4
#F = x^4 - 2*x^3*y + 14*x^2*y^2 - 16*x^2*y*z + 110*x^2*z^2 - 13*x*y^3 + 16*x*y^2*z - 110*x*y*z^2 + 52*y^4 - 4*y^3*z + 1199*y^2*z^2 + 3905*z^4
#F = 4*x^4 - 8*x^3*z + 92*x^2*y^2 - 148*x^2*y*z - 408*x^2*z^2 - 92*x*y^2*z + 148*x*y*z^2 + 412*x*z^3 - 371*y^4 - 742*y^3*z + 2991*y^2*z^2 + 302*y*z^3 - 2051*z^4
#a = 2
#b = 3
#F = b*z^4 + z*(x^3+y^3) + a*z^2*x*y + x^2*y^2
## The main functionality:
#X = Curve(P2.subscheme(F))

# The main functionality:
End = EndomorphismData(X, 100, have_oldenburg = True)
print End.period_matrix()
print End.geometric_representations()
print End.endomorphism_field()
print End.geometric().representations()
print End.over_base().representations()
K.<r> = NumberField(x^2 - 2)
print End.over_field(K).representations()
print End.over_field(K).structure()
print End.over_field(K).description()
print End.lattice()
print End.lattice().representations()
print End.lattice().structures()
print End.lattice().descriptions()
print End.verify_saturated()
print End.over_field(K).description()
print End.over_field(K).pretty_print()
print sagify_description(End.lattice().descriptions())
print End.lattice().pretty_print()

# Verification of decomposition data:
#Dec = End.decomposition()
#print Dec
#print Dec.field_of_definition()
#print Dec.impotents()
#print Dec.factors()
#print Dec.verify()

# Verification of geometric endomorphisms:
# Should have multiple parts
#print End.verify()

