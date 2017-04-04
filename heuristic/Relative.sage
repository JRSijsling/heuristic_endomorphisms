"""
 *  Relative splitting fields with GP optimizations
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
"""

def Relative_Splitting_Field(fs, bound = 0):
    bound_set = (bound != 0)
    F = magma.BaseRing(fs[1])
    overQQ = (magma.Degree(F) == 1)
    if overQQ:
        R.<x> = PolynomialRing(QQ)
    fs = sorted(fs, key = lambda f : -magma.Degree(f))
    K = F
    for f in fs:
        if not magma.HasRoot(f, K):
            for tup in magma.Factorization(f):
                K = magma.ExtendRelativeSplittingField(K, F, tup[1])
                if overQQ:
                    g = magma.DefiningPolynomial(K)
                    g = R(str(gp.polredabs(R(magma.Eltseq(g)))))
                    K = magma.NumberField(g)
                else:
                    K = magma.ClearFieldDenominator(K)
                if bound_set and magma.Degree(K) >= bound:
                    K = magma.DefineOrExtendInfinitePlaceFunction(K);
                    return K
    K = magma.DefineOrExtendInfinitePlaceFunction(K);
    return K
