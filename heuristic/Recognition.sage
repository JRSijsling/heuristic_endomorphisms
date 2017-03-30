"""
 *  Routines for recognizing complex numbers as elements of number fields
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
"""

def Common_Splitting_Field(AsPol, Bound = 0):
    # Input:   A collection of matrices AsPol with polynomial entries.
    # Output:  A common splitting field, represented by a tuple, in optimized
    #          form.
    # An optional bound makes this considerably faster to use in genus 2.
    R.<x> = PolynomialRing(QQ)
    Pols = reduce(lambda x,y: x + y, [A.sage().list() for A in AsPol])
    Pols = [ R(pol) for pol in Pols ]
    Pols_used = []
    f = x
    K = magma.RationalsAsNumberField()
    Bound_Set = (Bound != 0)
    # Sort by degree to find one likely to give a generator of the splitting
    # field:
    Pols = sorted(Pols, key = lambda pol : -pol.degree())
    for pol in Pols:
        if not magma.HasRoot(pol, K):
            Pols_used.append(pol)
            K = magma.MySplittingField(Pols_used)
            f = magma.DefiningPolynomial(K)
            cs = magma.Eltseq(f)
            f = R(str(gp.polredabs(R(cs))))
            K = magma.NumberField(f)
            if Bound_Set and magma.Degree(f) >= Bound:
                break
    return f.list()
