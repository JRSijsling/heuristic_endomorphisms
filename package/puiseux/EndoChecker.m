/***
 *  Endomorphism checking, Magma code
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
 */

function DiffBasis(X, is_hyp, is_plane_quartic);
/*
 * Input:   A curve X,
 *          whether or not X is hyperelliptic,
 *          and whether or not X is a plane quartic.
 * Output:  A basis of global differentials on X, represented by elements of
 *          the ambient.
 */

R := CoordinateRing(Ambient(X));
if is_hyp then
    x := R.1;
    y := R.2;
    f := DefiningEquations(X)[1];
    g := (Degree(f) - 1) div 2;
    c2 := MonomialCoefficient(f, y^2);
    c1 := &+[ MonomialCoefficient(f, x^i*y) * x^i : i in [0..g] ];
    K<x,y> := FunctionField(X);
    /* NOTE: these expressions are not particularly nice when called in Magma,
     * which leads to a slowdown. Still used for readability, and in fact this
     * is a responsibility of Magma to do well. */
    return [ x^(i-1) / (y + c1/(2*c2)) * Differential(x) : i in [1..g] ];
elif is_plane_quartic then
    x := R.1;
    y := R.2;
    f := DefiningEquations(X)[1];
    g := 3;
    K<x,y> := FunctionField(X);
    return [ (n / f) * Differential(x) : n in [x,y,1] ];
else
    // TODO: Provide names to variables in general case.
    K<x,y> := FunctionField(X);
    return BasisOfDifferentialsFirstKind(X);
end if;

end function;


function CandidateDivisors(X, g, is_hyp, is_planar, d);
/*
 * Input:   A curve X,
 *          its genus g,
 *          information about whether it is hyperelliptic and/or planar,
 *          and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

A := Ambient(X);
dimA := Dimension(A);
R := CoordinateRing(A);
F := BaseRing(R);
g := Genus(X);

if is_hyp then
    x,y := Explode(GeneratorsSequence(R));
    Rprod<x1, y1, x2, y2> := PolynomialRing(F, dimA * g);
    Xdivs := [ x^i : i in [0..(d div 2)] ] cat [ x^i*y : i in [0..((d - g - 1) div 2)] ];
elif is_planar then
    x,y := Explode(GeneratorsSequence(R));
    Rprod<x1, y1, x2, y2> := PolynomialRing(F, dimA * g);
    f := DefiningEquations(X)[1];
    Xdivs := [ x^i*y^j : i in [0..d], j in [0..(Degree(f) - 1)] ];
else
    // Despite the redundancy this usually seems to work well as a basis.
    // TODO: Name issues.
    Rprod<x1, y1, x2, y2> := PolynomialRing(F, dimA * g);
    Xdivs := &cat[ [ mon : mon in MonomialsOfDegree(R, n) ] : n in [1..d] ];
end if;

hs := [ hom<R -> Rprod | [ Rprod.j : j in [ ((i-1)*dimA + 1)..i*dimA ] ]> : i in [1..g] ];
CP := CartesianPower(Xdivs, g);
return [ &*[ hs[i](tup[i]) : i in [1..g] ] : tup in CP ];

end function;


function IrrCompsFromBranch(X, fs, n, P, Q);
/*
 * Input:   An embedded curve X,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Q.
 * Output:  The irreducible components corresponding to the kernel of N.
 */
// TODO: take more than one branch Q? Should probably use all of them in
// general. Later.

/* Recovering a linear system: */
e := Maximum([ Denominator(Valuation(c)) : c in Q ]);
M := [ ];
for f in fs do
    ev := Evaluate(f, P cat Q);
    // TODO: Check precision here.
    r := [ Coefficient(ev, i/e) : i in [0..n-1] ];
    Append(~M, r);
end for;
M := Matrix(M);
B := Basis(Kernel(M));

/* Coerce back to ground field (possible because of echelon form): */
F := BaseRing(X);
B := [ [ F ! c : c in Eltseq(b) ] : b in B ];

/* Corresponding equations: */
DEs := DefiningEquations(X);
R := Parent(DEs[1]);
Rprod := Parent(fs[1]);
d := #GeneratorsSequence(R);
g := #GeneratorsSequence(Rprod) div d;
hs := [ hom<R -> Rprod | [ Rprod.j : j in [ ((i-1)*d + 1)..i*d ] ]> : i in [1..g] ];
eqs := [ &+[ b[i] * fs[i] : i in [1..#fs] ] : b in B ];
eqs := eqs cat [ h(DE) : h in hs, DE in DEs ];

/* Corresponding scheme: */
A := AffineSpace(Rprod);
S := Scheme(A, eqs);

/* NOTE: The final two steps may be a time sink. Currently skipped. */
return [ S ];
Is := IrreducibleComponents(S);
return [ ReducedSubscheme(I) : I in Is ];

end function;


function IrrCompCheck(I, P0);
/*
 * Input:   An irreducible scheme I
 *          and a base point P0.
 * Output:  Whether or not I intersects P0 x X with the correct multiplicity at
 *          P0 and nowhere else.
 */

if Dimension(I) ne 1 then
    return false;
end if;
A := Ambient(I);
R := CoordinateRing(A);
g := #Eltseq(P0);
eqs := [ R.i - P0[i] : i in [ 1..g ] ];
S := Scheme(A, DefiningEquations(I) cat eqs);
test := (Dimension(S) eq 0) and (Degree(S) eq g) and (Degree(ReducedSubscheme(S)) eq 1);
return test;

end function;
        

function BasePointNonWeierstrassG2(X, As);

f, h := HyperellipticPolynomials(X);
n0 := 0;
while true do
    ev0 := Evaluate(f, n0);
    if ev0 ne 0 then
        break;
    end if;
    n0 +:= 1;
end while;
K := BaseRing(X);
R<t> := PolynomialRing(K);
L := SplittingField(t^2 - ev0);
P0 := [ L ! n0, Roots(t^2 - ev0, L)[1][1] ];
AsL := [ ChangeRing(A, L) : A in As ];
return P0, AsL;

end function;
