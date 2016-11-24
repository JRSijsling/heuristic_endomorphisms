intrinsic NonWeierstrassBasePointHyp(X::Crv, K::Fld, As::SeqEnum : B := 2^10) -> Crv, Pt, SeqEnum
{Given a hyperelliptic curve X, a field extension K, and a set of matrices As
over K, returns a base change of X, a point on that base change, and the base
changes of the matrices in As.}

/* We could look for points in the extension, but that is laborious */
Pts := RationalPoints(X : Bound := B);
Pts_nW := [ P : P in RationalPoints(X : Bound := B) | P[2]*P[3] ne 0 ];
if #Pts_nW ne 0 then
    Hts := [ Maximum([ Height(c) : c in Eltseq(P) ]) : P in Pts ];
    min, ind := Minimum(Hts);
    if Degree(K) eq 1 then
        K := Rationals();
    end if;
    XK := ChangeRing(X, K);
    PK := XK ! Pts_nW[ind];
    AsK := As;
    return XK, PK, AsK;
end if;

/* Find non-Weierstrass point: */
f := HyperellipticPolynomials(X);
n0 := 0;
while true do
    ev0 := Evaluate(f, n0);
    if ev0 ne 0 then
        break;
    end if;
    n0 +:= 1;
end while;

/* Extend and embed: */
R<t> := PolynomialRing(K);
L := SplittingField(t^2 - ev0);
if Degree(L) eq 1 then
    L := Rationals();
end if;
XL := ChangeRing(X, L);
PL := XL ! [ L ! n0, Roots(t^2 - ev0, L)[1][1] ];
AsL := [ ChangeRing(A, L) : A in As ];
if #Pts ne 0 then
    return XL, PL, AsL;
end if;

end intrinsic;

