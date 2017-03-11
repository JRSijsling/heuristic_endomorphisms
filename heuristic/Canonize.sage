"""
 *  Canonizing number fields
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

# Canonizes the number fields obtained in the endomorphism representation;
# the input and the output of these functions are all in string form, as by
# exporting from Magma to Sage we have already "broken" the link with actual
# objects.

# TODO: This is extremely ugly and needs an overhaul

def Minimal_Sequence(Ls):
    # Finds the sequence of minimal height in a list of sequences; if there are
    # two, then it return the first in the lexicographical ordering.
    sortf = lambda L : [ max([ QQ(c).height() for c in L ]), L[0] ]
    return sorted(Ls, key = sortf)

def Canonize_Field(frep):
    # Canonizes the number field associated to frep. Second return value gives
    # isomorphism, in the sense that it yields the roots of the canonized
    # polynomial in the original field.
    R.<x> = PolynomialRing(QQ)
    den = frep[len(frep) - 1]
    frepdiv = [ c / den for c in frep]
    f = R(frepdiv)
    fcan = R(str(gp.polredabs(f)))
    L.<r> = NumberField(f)
    froots = fcan.roots(L)
    return [ fcan.list(), [ rt[0].list() for rt in froots ] ]

def Canonize_Subfield(grep, frep):
    # Canonizes subfield given by grep in ambient frep that is supposed to have
    # been polredabs'ed already.
    R.<x> = PolynomialRing(QQ)
    f = R(frep)
    L.<r> = NumberField(f)
    # TODO: Sage sometimes bugs here when 1 is given, so we do that by hand:
    if L(grep) == L(1):
        return [ [0, 1], L(0).list() ]
    tup = L.subfield(L(grep))
    K = tup[0]
    phi = tup[1]
    can = Canonize_Field(K.0.minpoly().list())
    Kcan = can[0]
    rtreps = [ phi(K(rt)).list() for rt in can[1] ]
    return [ Kcan, Minimal_Sequence(rtreps)[0] ]

def Canonize_Subfield_And_Idempotents(grep, frep, idems):
    # Canonizes subfield given by grep in ambient frep that is supposed to have
    # been polredabs'ed already. Plus idempotents.
    R.<x> = PolynomialRing(QQ)
    f = R(frep)
    L.<r> = NumberField(f)
    # TODO: Sage sometimes bugs here when 1 is given, so we do that by hand:
    if L(grep) == L(1):
        return [ [0, 1], L(0).list(), magma.CanonizeMatrices(grep, frep, idems[1]) ]
    tup = L.subfield(L(grep))
    K = tup[0]
    phi = tup[1]
    can = Canonize_Field(K.0.minpoly().list())
    Kcan = can[0]
    rtreps = [ phi(K(rt)).list() for rt in can[1] ]
    rt = Minimal_Sequence(rtreps)[0]
    return [ Kcan, rt, magma.CanonizeMatrices(rt, frep, idems[1]) ]
