/***
 *  Verifying correspondences
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


// TODO: Make something similar for the planar case, Bound search only in case
// over QQ, correct field extension, allow for infty by patch, relative splitting field

intrinsic NonWeierstrassBasePointHyperelliptic(X::Crv, K::Fld : Bound := 2^10) -> SeqEnum
{Returns a non-Weierstrass point over a small extension of K.}

f, h := HyperellipticPolynomials(X);
g := 4*f + h^2; Y := HyperellipticCurve(g);

/* We could look for points in the extension, but that typically takes too much
 * time, so we restrict to the base */
Qs := RationalPoints(Y : Bound := Bound);
Qs_nW := [ Q : Q in Qs | not IsWeierstrassPlace(Place(Q)) ];
Qs_nW := [ Q : Q in Qs_nW | Q[3] ne 0 ];

if #Qs_nW ne 0 then
    Hts := [ Maximum([ Height(c) : c in Eltseq(Q) ]) : Q in Qs_nW ];
    min, ind := Minimum(Hts);
    L := K;
    Q := Qs_nW[ind];
    P := [ Q[1], (Q[2] - Evaluate(h, Q[1]))/2 ];

else
    /* Find non-Weierstrass point: */
    n0 := 0;
    while true do
        ev0 := Evaluate(g, n0);
        if ev0 ne 0 then
            break;
        end if;
        n0 +:= 1;
    end while;

    /* Extend and embed: */
    R<t> := PolynomialRing(K);
    L := SplittingField(t^2 - ev0);
    Q := [ L ! n0, Roots(t^2 - ev0, L)[1][1] ];
    P := [ Q[1], (Q[2] - Evaluate(h, Q[1]))/2 ];
end if;

return P;

end intrinsic;


intrinsic NonWeierstrassBasePointPlane(X::Crv, K::Fld : Bound := 2^10) -> SeqEnum
{Returns a non-Weierstrass point over a small extension of K.}

// TODO: We can do a bit better in the case where there is a rational point by
// taking lines through it.

Ps := RationalPoints(X);
Ps_nW := [ P : P in Ps | not IsWeierstrassPlace(Place(P)) ];
if #Ps_nW ne 0 then
    Hts := [ Maximum([ Height(c) : c in Eltseq(P) ]) : P in Ps_nW ];
    min, ind := Minimum(Hts);
    L := K;
    P := Ps_nW[ind];
    return P;
end if;

f := DefiningPolynomial(X);
R<x,y,z> := PolynomialRing(K, 3);
S<t> := PolynomialRing(K);

/* If there is a rational point, then take lines through it: */
if #Ps ne 0 then
    P0 := Ps[1];
    x0, y0, z0 := Explode(Eltseq(P0));
    n0 := 0;
    while true do
        if z0 ne 0 then
            h := hom< R -> S | [ n0*t + x0, t + y0, z0 ]>;
        else
            h := hom< R -> S | [ n0*t + x0, y0, t + z0 ]>;
        end if;
        Fac := Factorization(h(f));
        for tup in Fac do
            fac := tup[1];
            L := NumberField(fac);
            rt := Roots(fac, L)[1][1];
            XL := ChangeRing(X, L);
            if z0 ne 0 then
                P := XL ! [ n0*rt + x0, rt + y0, z0 ];
            else
                P := XL ! [ n0*rt + x0, y0, rt + z0 ];
            end if;
            if not IsWeierstrassPlace(Place(P)) then
                return P;
            end if;
        end for;
        n0 +:= 1;
    end while;
end if;

/* Otherwise lines through (0,0,1): */
n0 := 0;
while true do
    h := hom< R -> S | [ n0, t, 1 ]>;
    Fac := Factorization(h(f));
    for tup in Fac do
        L := NumberField(tup[1]);
        rt := Roots(h(f), L)[1][1];
        XL := ChangeRing(X, L);
        P := XL ! [ n0, rt, 1 ];
        if not IsWeierstrassPlace(Place(P)) then
            return P;
        end if;
    end for;
    n0 +:= 1;
end while;

end intrinsic;


intrinsic NonWeierstrassBasePoint(X::Crv, K::Fld : Bound := 2^10) -> SeqEnum
{Returns a non-Weierstrass point over a small extension of K.}

if Type(X) eq CrvHyp then
    return NonWeierstrassBasePointHyperelliptic(X, K : Bound := Bound);
elif Type(X) eq CrvPln then
    return NonWeierstrassBasePointPlane(X, K : Bound := Bound);
else
    error "Not implemented for general curves yet";
end if;

end intrinsic;


intrinsic Correspondence(X::Crv, P::Pt, Y::Crv, Q::Pt, A::.) -> .
{Gives certificate.}
// TODO: Update this after changing DivisorFromMatrixSplit
// TODO: Add decent compositum (currently we assume that Q is defined over the field for P)
// TODO: Carry P around or make it part fof object?
// TODO: Correct matrices (we have transformed the curves). This should arguably happen before.

L := Parent(P[1]);
XL := ChangeRing(X, L); YL := ChangeRing(Y, L);
PL := XL ! P; QL := YL ! Q;
AL := ChangeRing(A, L);

if IsScalar(AL) then
    return true, "Scalar: OK";
else
    // TODO: Add theoretical lower bounds (but for now we can experiment)
    LowerBound := 1;
    return true, DivisorFromMatrixSplit(XL, PL, YL, QL, Transpose(AL) : LowerBound := LowerBound);
end if;

end intrinsic;
