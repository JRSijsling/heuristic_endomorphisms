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

Xs = [ ]

# RM over QQ
f = x^8 + 14*x^4 + 1
h = R(0)
Xs.append(HyperellipticCurve(f, h))

for X in Xs:
    print X
    # The main functionality
    Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

    #print "Period matrix:"
    #print Endo.period_matrix()

    print "Field of definition:"
    print Endo.endomorphism_field()

    #print "Testing Rosati and degree bound:"
    #A = Endo._geo_rep_list_[1][1]
    #print A
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

EndoGeo = Endo.geometric()
EndoDict = EndoGeo.full()
alg = EndoDict['algebra']
A = alg['alg_QQ']
print A
