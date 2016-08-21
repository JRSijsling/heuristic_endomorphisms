"""
 *  Debugging functionality over a finite field
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

load('EndoChecker.sage')
load('LocalInfo.sage')
magma.load('EndoChecker.m')
magma.load('LocalInfo.m')

# TODO: Use split primes only for speedup?
ps = [ n for n in [10^6..10^6 + 100] if n.is_prime() ]
for p in ps:
    F = FiniteField(p)
    R.<x> = PolynomialRing(F)
    
    f = 1 + 2*x + 7*x^2 + 4*x^3 + 7*x^4 + 2*x^5 + x^6
    X = HyperellipticCurve(f)
    P0 = [0, 1]
    M = Matrix(F, [[3, 0], [0, 3]])
    M = Matrix(F, [[0, 1], [1, 0]])
    M = Matrix(F, [[0, -1], [-1, 0]])

    EC = EndoChecker(X, P0)
    n = 300
    print p
    print EC.LI.nth_approxs(M, n)
    print EC.LI.nth_approxs_old(M, n)
