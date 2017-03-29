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
P2.<x,y,z> = ProjectiveSpace(QQ, 2)
F = x^4 - 5*x^3*y + 3*x^3*z + 3*x^2*y^2 + 4*x^2*y*z + x^2*z^2 + 2*x*y^3 + 4*x*y^2*z + x*y*z^2 - x*z^3 - 2*y^4 + 4*y^3*z - 2*y^2*z^2 + 3*y*z^3 - 2*z^4

# The main functionality:
X = Curve(P2.subscheme(F))
End = EndomorphismData(X, prec = 300)
print End.period_matrix()
print End.geometric_representations()
print End.field_of_definition()
print End.lattice()
#print End.geometric().description()
#print End.over_base().description()

# Verification of decomposition data:
#Dec = End.decomposition()
#print Dec
#print Dec.field_of_definition()
#print Dec._idems_[1]
#print Dec.factors()
#print Dec.certificate_g2()

# Verification of geometric endomorphisms:
#print End.geometric_representations_check()

