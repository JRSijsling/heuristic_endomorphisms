"""
 *  Endomorphism checking, Sage code
 *
 *  Copyright (C) 2016  J.R. Sijsling (sijsling@gmail.com)
 *
 *  Distributed under the terms of the GNU General License (GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc., 51
 *  Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

class EndoChecker:
    def __init__(self, X, P0, degree_bound = 0):
        X = magma(X)
        P0 = X(magma(P0))
        self.degree_bound = degree_bound
        self.is_hyp = magma.IsHyperellipticCurve(X)
        self.is_planar = magma.IsPlaneCurve(X)
        self.is_smooth = magma.IsSingular(X)
        self.is_plane_quartic = (self.is_smooth) and (magma.Genus(X) == 3)
        self.F = magma.BaseRing(X)
        # Taking affine patch to get divisors later:
        self.X, self.P0 = magma.AffinePatch(X, P0, nvals = 2)
        assert not magma.IsWeierstrassPlace(magma.Place(self.P0))
        self.A = magma.Ambient(self.X)
        self.S = magma.CoordinateRing(self.A)
        self.B = magma.DiffBasis(self.X, self.is_hyp, self.is_plane_quartic)
        self.g = len(self.B)
        self.LI = LocalInfo(self.X, self.P0, self.is_planar, self.B)

    def __repr__(self):
        return "Integration of endomorphisms from tangent matrices for the {}. Base point used: {}.".format(str(self.X), str(self.P0))

    def candidate_divisors(self, d):
        return magma.CandidateDivisors(self.X, self.g, self.is_hyp, self.is_planar, d)

    def divisor_from_matrix(self, M, margin = 2^4):
        # We start at a suspected estimate and then increase degree until we find an appropriate divisor:
        d = self.degree_bound
        while True:
            print "Trying degree {}...".format(str(d))
            fs = self.candidate_divisors(d)
            n = len(fs) + margin
            print "Precision required: {}.".format(str(n))
            # Take non-zero image branch:
            print "Lifting the point..."
            P, alphaP = self.LI.nth_approxs(M, n)
            print "done lifting."
            for Q in alphaP:
                if not all([ c == 0 for c in Q ]):
                    break
            # Fit a divisor to it:
            ICs = magma.IrrCompsFromBranch(self.X, fs, n, P, Q)
            for I in ICs:
                if magma.IrrCompCheck(I, self.P0):
                    print "Divisor found in degree {}!".format(str(d))
                    return I
            # If that does not work, give up and try one degree higher:
            d += 1
