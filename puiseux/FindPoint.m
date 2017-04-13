/***
 *  Finding a suitable point on a curve
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

// TODO: Make something similar for the planar case, Bound search only in case
// over QQ, correct field extension, allow for infty by patch

intrinsic NonWeierstrassBasePointHyperelliptic(X::Crv, As::SeqEnum : Bound := 2^10) -> Crv, SeqEnum, Pt
{Given a hyperelliptic curve X, and a set of matrices As over an extension of
the base field of X, returns a base change of X, a point on that base change,
and the base changes of the matrices in As.}

K := BaseRing(As[1]);
f, h := HyperellipticPolynomials(X);
g := 4*f + h^2; Y := HyperellipticCurve(g);

/* We could look for points in the extension, but that typically takes too much
 * time, so we restrict to the base */
Qs := RationalPoints(Y : Bound := Bound);
Qs_nW := [ Q : Q in RationalPoints(Y : Bound := Bound) | not IsWeierstrassPlace(Place(Q)) ];
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

XL := ChangeRing(X, L); YL := ChangeRing(Y, L);
PL := XL ! P; QL := YL ! Q;
AsL := [ ChangeRing(A, L) : A in As ];
//return XL, AsL, PL;
return YL, AsL, QL;

end intrinsic;


intrinsic NonWeierstrassBasePointPlane(X::Crv, As::SeqEnum : Bound := 2^10) -> Crv, Pt, SeqEnum
{Given a curve X, a field extension K, and a set of matrices As over K, returns
a base change of X, a point on that base change, and the base changes of the
matrices in As.}

return 0;

end intrinsic;


intrinsic NonWeierstrassBasePoint(X::Crv, As::SeqEnum : Bound := 2^10) -> Crv, Pt, SeqEnum
{Given a curve X, a field extension K, and a set of matrices As over K, returns
a base change of X, a point on that base change, and the base changes of the
matrices in As.}

if Type(X) eq CrvHyp then
    return NonWeierstrassBasePointHyperelliptic(X, As : Bound := Bound);
elif Type(X) eq CrvPln then
    return NonWeierstrassBasePointPlane(X, As : Bound := Bound);
end if;

end intrinsic;
