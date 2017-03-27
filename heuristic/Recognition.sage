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
            #print "Calculating splitting field of:"
            #print Pols_used
            K = magma.MySplittingField(Pols_used)
            f = magma.DefiningPolynomial(K)
            cs = magma.Eltseq(f)
            # Apply GP simplification plus a dirty hack because I do not know
            # better:
            # NOTE: Direct conversion of f fails completely!
            f = R(str(gp.polredabs(R(cs))))
            K = magma.NumberField(f)
            if Bound_Set and magma.Degree(f) >= Bound:
                break
    return f.list()


# What follows is no longer used; the more hands-on Magma implementation that
# uses LLL seems better.

def Polynomialize_Element(a, epscomp = epscomp):
    # Input:    An element of a complex field.
    # Output:   A minimal polynomial one of whose roots approximates a well.
    Ca = a.parent()
    d = 1
    h = h0 = 2^4
    # Note that the variable of the parent is important; we cannot coerce algdep
    # otherwise.
    R.<x> = PolynomialRing(ZZ)
    while true:
        p = algdep(a, d, height_bound = h)
        if p in R:
            for tup in p.factor():
                q = tup[0]
                if q(a).abs() < epscomp:
                    return q
        d += 1
        h *= h0
        # Making sure of loop exit:
        if d > 100:
            raise ValueError("No polynomials of small degree found")

def Polynomialize_Matrix(A, epscomp = epscomp):
    # Input:    A matrix over a complex field.
    # Output:   The same matrix with its entries replaced by the polynomials
    #           obtained by running the previous algorithm.
    return Matrix([[Polynomialize_Element(x, epscomp = epscomp) for x in r] for
        r in A.rows()])
