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

load("../Initialize.sage")

# Hyperelliptic tests over QQ
R.<x> = PolynomialRing(QQ)

# Curve input: specify g and h in its equation y^2 + h y = g.

# RM over QQ
f = -x^5
h = x^3 + x + 1
# Stupid split case
f = x^6 + x^2 + 1
h = R(0)
# Genus 2 factor
f = x^8 - 12*x^7 + 50*x^6 - 108*x^5 + 131*x^4 - 76*x^3 - 10*x^2 + 44*x - 19
h = R(0)
# Splitting functionality
f = 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x
h = x
# FindPoint
f = 3*x^3 - 2*x^2 + 6*x + 2
h = x^3 + x
# Subfield
f = x^4 + x^3 + 2*x^2 + x + 1
h = x^3 + x^2 + x + 1
# Wrong point
f = x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8
h = 0
# Largest genus 2 case
f = x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1
h = R(0)
# Test case
f = x^5 + x^4 + 2*x^3 + x^2 + x
h = x^2 + x
# QM
f = x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1
h = R(0)

# Hyperelliptic tests over extension
#F.<r> = NumberField(x^2 - 5)
#R.<x> = PolynomialRing(F)
#f = x^5 + r*x^3 + x
#h = R(0)

X = HyperellipticCurve(f, h)

## Plane tests
#P2.<x,y,z> = ProjectiveSpace(QQ, 2)
#F = x^4 - 3*x^3*y + 3*x^3*z + 3*x^2*y^2 + 4*x^2*y*z + x^2*z^2 + 2*x*y^3 + 4*x*y^2*z + x*y*z^2 - x*z^3 - 2*y^4 + 4*y^3*z - 2*y^2*z^2 + 3*y*^3 - 2*z^4
#F = x^4 - 2*x^3*y + 14*x^2*y^2 - 16*x^2*y*z + 110*x^2*z^2 - 13*x*y^3 + 16*x*y^2*z - 110*x*y*z^2 + 52*y^4 - 4*y^3*z + 1199*y^2*z^2 + 3905*z^4
#F = 4*x^4 - 8*x^3*z + 92*x^2*y^2 - 148*x^2*y*z - 408*x^2*z^2 - 92*x*y^2*z + 148*x*y*z^2 + 412*x*z^3 - 371*y^4 - 742*y^3*z + 2991*y^2*z^2 + 302*y*z^3 - 2051*z^4
#a = 2
#b = 3
#F = b*z^4 + z*(x^3+y^3) + a*z^2*x*y + x^2*y^2
## The main functionality
#X = Curve(P2.subscheme(F))

# LMFDB tests:
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -1, 0, 0, 0, 1]), R([]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([-2, -8, -10, -5, 0, 1]), R([0, 0, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([-2]), R([0, 0, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 0, 0, 0, 0, 1]), R([1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 0, 0, 0, 0, -1]), R([1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([1, 3, 6, 7, 6, 3, 1]), R([0, 1, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -1, 0, -1]), R([1, 1, 1, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 1, 2, 2, 1]), R([1, 1, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -1, 1, 1, -3, 2]), R([1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 0, 0, 1, 1]), R([1, 1, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([-1, 0, 2, 0, -2, 0, 1]), R([]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -1, 0, -1, 0, 1]), R([]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([-1]), R([1, 0, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -1, 1, 0, -1, 1]), R([]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([1, 0, 0, 1]), R([1, 0, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 0, 0, 0, 1]), R([1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([-7, 0, 0, 0, 1]), R([0, 1, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 1, 0, 1]), R([1, 0, 0, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([1, 1, 2, 1, 1]), R([1, 1, 1, 1]))
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 1, 1]), R([1, 0, 0, 1]))
X = C

# TODO: Algebraization problem
#R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -4, 6, 0, -3, 1]), R([]))

# The main functionality
Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

#print "Period matrix:"
#print Endo.period_matrix()

print "Field of definition:"
print Endo.endomorphism_field()

#print "Testing Rosati and degree bound:"
#A = Endo._geo_rep_list_[2][1]
#print Endo.rosati_involution(A)
#print Endo.degree_estimate(A)

print "Over several fields:"
#print Endo.geometric().representation()
#print Endo.over_base().representation()
K.<s> = NumberField(x^2 - 2)
overK = Endo.over_field(K)
print K
#print overK.representation()
#print overK.algebra()
#print overK.description()
print overK.pretty_print()

print "Examples of lattices:"
#print Endo.lattice()
#print Endo.lattice().representations()
#print Endo.lattice().algebras()
#print Endo.lattice().descriptions()
print Endo.lattice().pretty_print()

#print "Verification:"
#print Endo.verify()

#print "Decomposition:"
#Dec = Endo.decomposition()
#print Dec
#print Dec.decomposition_field()
#print Dec.idempotents()
#print Dec.projections()
#print Dec.factors()
#print Dec.verify()
#print Dec.morphisms()
