"""
 *  Determining an elliptic curve from a lattice
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

def Elliptic_Curve_From_Lattice(L, frep, h, prec = prec, epscomp = epscomp,
        epsLLL = epsLLL): 
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

# Part of TODO:
# This allows us to find the curves. To get points later, use
# gp.ellwp(w, 0.00001), which indeed has a quadratic pole.
