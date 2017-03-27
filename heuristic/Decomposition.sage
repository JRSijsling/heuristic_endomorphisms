"""
 *  Determining an elliptic curve from a lattice
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
"""

def Elliptic_Curve_From_Lattice(L, frep, h, prec = prec, epscomp = epscomp, epsLLL = epsLLL):
    # Input:    A lattice of complex numbers L, a representation of a number
    #           field frep, and number representing a homomorphism h from that
    #           field into the complex numbers.
    # Output:   An algebraic equation of the corresponding elliptic curve.
    gp.set_precision(prec)
    g4 = gp.elleisnum(L, 4, flag = 1)
    g6 = gp.elleisnum(L, 6, flag = 1)
    eisvals_an = [ 12 * g4 , 216 * g6 ]
    eisvals_an = [ eisval.sage() for eisval in eisvals_an ]
    # TODO: Next line can be replaced by a smaller field; if we just keep track
    # of its hom then that gives descents.
    eisvals_alg = [ magma.AlgebraizeElementInField(magma(eisval), magma(frep),
        magma(h), epscomp = epscomp, epsLLL = epsLLL) for eisval in eisvals_an ]
    return eisvals_alg
