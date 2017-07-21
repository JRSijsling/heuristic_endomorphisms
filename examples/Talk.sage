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
F = QQ
R.<x> = PolynomialRing(F)
Xs = [ ]

# RM over QQ
f = -x^5
h = x^3 + x + 1
Xs.append(mHyperellipticCurve(f, h))

# Splits over extension
f = x^4 + x^3 + 2*x^2 + x + 1
h = x^3 + x^2 + x + 1
Xs.append(mHyperellipticCurve(f, h))

# QM
f = x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1
h = R(0)
Xs.append(mHyperellipticCurve(f, h))


# Hyperelliptic test over extension
R.<t> = PolynomialRing(QQ)
F.<r> = NumberField(t^2 - 5)
R.<x> = PolynomialRing(F)

f = x^5 + r*x^3 + x
h = R(0)
Xs.append(mHyperellipticCurve(f, h))


# Plane test
F = QQ
P2.<x,y,z> = ProjectiveSpace(F, 2)

f = y^3*z - (x^4 + 2*x^2*z^2 + 5*z^4)
Xs.append(mPlaneCurve(f))

f = x^4 - x^3*y + 2*x^3*z + 2*x^2*y*z + 2*x^2*z^2 - 2*x*y^2*z + \
    4*x*y*z^2 - y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4;
Xs.append(mPlaneCurve(f))


# Here we go
for X in Xs:
    print X
    print ""
    Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

    print "Field of definition:"
    print Endo.endomorphism_field()
    print ""

    print "Geometric representation:"
    print Endo.geometric().representation()
    print ""

    print "Lattice:"
    print Endo.lattice().pretty_print()
    print ""

    print "Decomposition:"
    Dec = Endo.decomposition()
    print Dec.idempotents()
    print Dec.factors()
    #print Dec.verify()
    print ""
