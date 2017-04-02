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
F = x^4 - 3*x^3*y + 3*x^3*z + 3*x^2*y^2 + 4*x^2*y*z + x^2*z^2 + 2*x*y^3 + 4*x*y^2*z + x*y*z^2 - x*z^3 - 2*y^4 + 4*y^3*z - 2*y^2*z^2 + 3*y*z^3 - 2*z^4
F = x^4 - 2*x^3*y + 14*x^2*y^2 - 16*x^2*y*z + 110*x^2*z^2 - 13*x*y^3 + 16*x*y^2*z - 110*x*y*z^2 + 52*y^4 - 4*y^3*z + 1199*y^2*z^2 + 3905*z^4
#F = 4*x^4 - 8*x^3*z + 92*x^2*y^2 - 148*x^2*y*z - 408*x^2*z^2 - 92*x*y^2*z + 148*x*y*z^2 + 412*x*z^3 - 371*y^4 - 742*y^3*z + 2991*y^2*z^2 + 3302*y*z^3 - 2051*z^4
a = 2
b = 3
F = b*z^4 + z*(x^3+y^3) + a*z^2*x*y + x^2*y^2

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

