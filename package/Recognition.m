/***
 *  Functionality for recognizing complex numbers as algebraic numbers and
 *  related optimizations
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

function PolynomializeElement(a : epscomp := epscomp0, epsLLL := epscomp0);
// Input:    An element of a complex field.
// Output:   A minimal polynomial one of whose roots approximates a well.
R<x> := PolynomialRing(Rationals());
Ca := Parent(a);
Ra := RealField(Precision(Ca));
// d is redundant but useful for clarity
d := 0;
h := 1;
h0 := 10;
powers := [ Ca ! 1 ];
while true do
    // Increase height:
    d +:= 1;
    h *:= h0;
    Append(~powers, powers[#powers] * a);
    M := Transpose(Matrix([powers]));
    // Now split and take an IntegralKernel:
    MSplit := SplitPeriodMatrix(M);
    Ker := IntegralKernel(MSplit : epsLLL := epsLLL);
    for row in Rows(Ker) do
        // Height condition:
        if &and[ Abs(c) le h : c in Eltseq(row) ] then
            pa := &+ [ row[i + 1] * x^i : i in [0..d] ];
            // Factor (to eliminate redundancy) and check:
            Fac := Factorization(pa);
            for tup in Fac do
                fa := tup[1];
                if Abs(Evaluate(fa, a)) lt epscomp then
                    return fa;
                end if;
            end for;
        end if;
    end for;
    // Making sure of loop exit:
    if d gt 100 then
        error Error("No polynomials of small degree found");
    end if;
end while;

end function;


function PolynomializeMatrix(A : epscomp := epscomp0, epsLLL := epscomp0);
// Input:    A matrix over a complex field.
// Output:   The same matrix with its entries replaced by the polynomials
//           obtained by running the previous algorithm.
    return Matrix([[ PolynomializeElement(a : epscomp := epscomp, epsLLL :=
        epsLLL) : a in Eltseq(r) ] : r in Rows(A) ]);
end function;


// The following function is generalized by AlgebraizeElementInField:
function FractionalApproximation(a : epscomp := epscomp0, epsLLL := epsLLL0);
    // Input:    An element of a complex field.
    // Output:   A fraction that approximates to the given precision
    //           (using continued fractions is likely much better).
    M := Matrix([[1], [-Real(a)]]);
    // A loop should follow, but not for now:
    K := IntegralKernel(M : epsLLL := epsLLL);
    q := K[1,1] / K[1,2];
    if Abs(q - a) lt epscomp then
        return q;
    else
        error Error("LLL not does return a sufficiently good fraction");
    end if;
end function;


function AlgebraizeElementInField(a, frep, h : epscomp := epscomp0,
    epsLLL := epsLLL0);
// Input:    A complex number a, a representation of a number field frep, and a
//           number representing a homomorphism h from that field into the
//           complex numbers.
// Output:   A power basis representation of x as an element of the field
//           corresponding to frep.
Ca := Parent(a);
ran := Ca!h;
Ra := RealField(Precision(Ca));
K := NumberField(Polynomial(frep));
d := Degree(K);
// FIXME: Usual hateful dichotomy. Maybe take roots of DefiningPolynomial or
// defining polynomial?
if d ne 1 then
    r := K.1;
else
    r := K ! 0;
end if;
// Creating algebraic and analytic powers:
// FIXME: Relies on 0^0 = 1.
powers := [ r^n : n in [0..d-1] ];
powersan := [ ran^n : n in [0..d-1] ];
// Column matrix of embedding of basis elements plus (minus) a:
Man := VerticalJoin(Transpose(Matrix([powersan])), Matrix([[-a]]));
// Now split and take an IntegralKernel:
ManSplit := SplitPeriodMatrix(Man);
Ker := IntegralKernel(ManSplit : epsLLL := epsLLL);
for R in Rows(Ker) do
    den := R[d + 1];
    if den ne 0 then
        s := (&+[ R[i] * powers[i] : i in [1..d] ])/ den;
        san := (&+[ R[i] * powersan[i] : i in [1..d] ]) / den;
        // Check correct to given precision:
        if Abs(san - a) lt epscomp then
            return Eltseq(s);
        end if;
    end if;
end for;

error Error("Fail to algebraize element in ambient");

end function;


// A slight variation of the above (in the application below, we could even
// replace it with AlgebraizeElementInField):
function NearbyRoot(a, f, h : epscomp := epscomp0);
// Input:   A polynomial f, a complex number a that is an approximation of a
//          complex root of f, and a homomorphism h of a number field K into the
//          parent of a.
// Output:  A root of f in K whose embedding under h is close to a.

K := Domain(h);
rs := [ tup[1] : tup in Roots(f, K) ];
for r in rs do
    if Abs(h(r) - a) lt epscomp then
        return r;
    end if;
end for;

error Error("Failed to find a nearby root");

end function;


function AlgebraizeMatricesInField(As, AsPol, frep : epscomp := epscomp0);
// Input:   A collection of matrices As and their polynomizations AsPol, along
//          with a number field, specified by a tuple frep.
// Output:  Algebraic representations of these matrices and representations of
//          the corresponding number field and embedding.

// Choose complex embedding; need to hack this because of Magma's stupidity.

f := Polynomial(frep);
prec := Precision(BaseRing(As[1]));
if Degree(f) eq 1 then
    // FIXME: Get rid of RationalsAsNumberField():
    K<r> := RationalsAsNumberField();
    h := hom<K -> ComplexField(prec) | 1>;
    fhom := h(0);
else
    K<r> := NumberField(f);
    ra := Roots(f, ComplexField(prec))[1][1];
    h := hom<K -> ComplexField(prec) | ra>;
    fhom := h(r);
end if;

// Unfortunately we cannot zip matrices in Magma, so we resort to running over
// indices:
L := #As;
M := #Rows(As[1]);
N := #Rows(Transpose(As[1]));
// TODO: This step may be too expensive, and we should use
// AlgebraizeElementInField
AsAlg := [Matrix(K, [[ NearbyRoot(As[l][m][n], AsPol[l][m][n], h : epscomp :=
    epscomp) : n in [1..N]] : m in [1..M] ]) : l in [1..L]];

// Return representation along with root needed to reconstruct:
return AsAlg, fhom;

end function;


function IntegralRepresentationNF(K);
// Input:   A number field K.
// Output:  A number field L with a small integral defining polynomial and an
//          isomorphism h from K to L.

d := Degree(K);
if d eq 1 then
    // In the case of degree 1 we use the rationals as a number field for
    // reasons of uniformity:
    L := RationalsAsNumberField();
    // FIXME: We skip the above:
    L := Rationals();
    if K eq Rationals() then
        return L, hom<K -> L | >;
    else
        return L, hom<K -> L | L ! (K.1)>;
    end if;
    L := Rationals();
else
    // We cannot quite apply OptimizedRepresentation straight away since it does
    // not always make the defining polynomial integral. Instead of this we find
    // a small element of the maximal order that defines the field and apply the
    // function to that instead:
    ZK := Integers(K);
    r := ZK.1;
    B := 0;
    // This loop is guaranteed to terminate:
    while true do
        g := MinimalPolynomial(r);
        if Degree(g) eq d then
            L := OptimizedRepresentation(NumberField(g));
            test, h := IsIsomorphic(K, L);
            return L, h;
        end if;
        r := ZK!([Random([-B..B]) : n in [1..d]]);
        B +:= 1;
    end while;
end if;

end function;


// Deprecated:
function PartialLMFDBLabel(K);
// Input:   A number field.
// Output:  The first three entries of its LMFDB label.

d := Degree(K);
if Degree(K) eq 1 then
    return [1, 1, 1];
else
    r := #[ v : v in InfinitePlaces(K) | IsReal(v) ];
    D := Abs(Discriminant(Integers(K)));
    return [d, r, D];
end if;

end function;
