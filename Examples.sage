"""
 *  Some examples of endomorphism rings of genus 2 hyperelliptic curves
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
"""

# Add if no initilization script set:
load('Initialize.sage')

# Ambient ring:
R.<x> = PolynomialRing(QQ)
#F.<r> = NumberField(x^2 - 2)
#R.<x> = PolynomialRing(F)

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
#f = 3*x^3 - 2*x^2 + 6*x + 2
#h = x^3 + x
# Debug (TODO: Subfield):
#f = x^4 + x^3 + 2*x^2 + x + 1
#h = x^3 + x^2 + x + 1
# Debug (TODO: Point returned is wrong):
#f = x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8
#h = 0
# Case where [RR, CC] occurs:
f = x^5 + x^4 + 2*x^3 + x^2 + x
h = x^2 + x
# Genus 3:
f = x^7 + x^6 + x^5 + x^3 + x^2 + x
h = x^4 + x^2 + 1

# The main functionality:
X = HyperellipticCurve(f, h)
End = EndomorphismData(X, 300, have_oldenburg = True)
print End.period_matrix()
print End.geometric_representations()
#print End.endomorphism_field()
#print End.lattice()
#print End.geometric().description()
#print End.over_base().description()
#print End.field_of_definition()

# Verification of decomposition data:
#Dec = End.decomposition()
#print Dec
#print Dec.field_of_definition()
#print Dec._idems_[1]
#print Dec.factors()
#print Dec.certificate_g2()

# Verification of geometric endomorphisms:
#print End.geometric_representations_check()

