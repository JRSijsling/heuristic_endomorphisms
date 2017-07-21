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
f = x^4 + x^3 + 3*x^2 + x + 2
h = x^3 + x^2 + x
f = 4*f + h^2
h = 0
Xs.append(mHyperellipticCurve(f, h))

for X in Xs:
    print X
    Endo = EndomorphismData(X, prec = 300, have_oldenburg = False)

    print "Representation on homology:"
    gens = Endo.geometric().full()['representation']
    homs = magma([ gen['homology'] for gen in gens ])
    R = homs[1] + homs[2]
    print R^2
