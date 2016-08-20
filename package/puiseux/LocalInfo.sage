"""
 *  Local information around a point, Sage code
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

class LocalInfo:
    def __init__(self, X, P0, is_planar, B):
        self.X = X
        self.P0 = P0
        self.B = B
        self.g = len(self.B)
        self.is_planar = is_planar
        # For normalizing the differential representation:
        self.Bnorm, self.T = magma.NormalizeDiffBasis(self.X, self.P0, self.is_planar, self.B, nvals = 2)

    def __repr__(self):
        return "Puiseux lifting on the {}. Base point used: {}. Basis of differentials: {}. Normalized basis of differentials: {}".format(str(self.X), str(self.P0), str(self.B), str(self.Bnorm))

    def nth_approxs(self, M, n):
        # Conjugate in order to work with the normalized basis:
        M = self.T * magma(M) * (self.T)^(-1)
        # NOTE: In the non-planar case the following functionality assumes that
        # Magma picks the same uniformizer throughout as I see no way to pass
        # it on:
        P = magma.DevelopInUnif(self.X, self.P0, self.is_planar, n + self.g)
        return [ P, magma.NthApproxs(self.X, P, self.is_planar, self.Bnorm, M, n) ]

    def nth_approxs_old(self, M, n):
        # Conjugate in order to work with the normalized basis:
        M = self.T * magma(M) * (self.T)^(-1)
        # NOTE: In the non-planar case the following functionality assumes that
        # Magma picks the same uniformizer throughout as I see no way to pass
        # it on:
        P = magma.DevelopInUnif(self.X, self.P0, self.is_planar, n + self.g)
        return [ P, magma.NthApproxsOld(self.X, P, self.is_planar, self.Bnorm, M, n) ]
