/***
 *  Factors from lattices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic FactorReconstructG1(P::., K::Fld) -> .
{Reconstructs elliptic curve factor from a period lattice.}

P := Eltseq(P);
CC := Parent(P[1]); RR := RealField(CC);
g4CC := 120 * (1/P[1])^4 * ZetaFunction(RR, 4) * Eisenstein(4, P);
g6CC := 280 * (1/P[1])^6 * ZetaFunction(RR, 6) * Eisenstein(6, P);
g4 := AlgebraizeElementInRelativeField(g4CC, K);
g6 := AlgebraizeElementInRelativeField(g6CC, K);
R<x> := PolynomialRing(K);
// TODO: Check that this is the correct curve (depends on differentials used)
return HyperellipticCurve([ 4*x^3 - g4*x - g6, R ! 0 ]);

end intrinsic;


intrinsic TwistDifferentialBasis(P::., X::.) -> .
{Twists the curve and differential basis as needed.}

K := BaseRing(X);
CC := Parent(P[1,1]); RCC := PolynomialRing(CC);
f, h := HyperellipticPolynomials(X);
// TODO: Move iota to keyword arguments?
fCC := EmbedAtInfinitePlace(f, K`iota, RCC); hCC := EmbedAtInfinitePlace(h, K`iota, RCC);
Q := PeriodMatrix(fCC, hCC);
AsRs := GeometricIsogenyBasisApproximations(P, Q); A := AsRs[1][1];

// TODO: In general we will need a function to test for lattice isomorphism. It
// will also need to be verified... also, we could get a triple extension here.
// The code is relatively ugly in order to avoid this. We would like to take
// another relative splitting field. Anyway.
pols := RelativeMinimalPolynomials(Eltseq(A), K);
sqrtCC := CC ! 1; d := 1;
for pol in pols do
    if Degree(pol) eq 2 then
        polCC := EmbedAtInfinitePlace(pol, K`iota, RCC);
        a := Coefficient(polCC, 2); b := Coefficient(polCC, 1); c := Coefficient(polCC, 0);
        aCC := Coefficient(polCC, 2); bCC := Coefficient(polCC, 1); cCC := Coefficient(polCC, 0);
        sqrtCC := Sqrt(bCC^2 - 4*aCC*cCC); d := b^2 - 4*a*c;
    end if;
end for;

WA := AlgebraizeMatrixInRelativeField(sqrtCC * A, K);
g := 4*f + h^2;
R<x> := Parent(g);
g := Evaluate(g, (WA[1,1]*x + WA[1,2])/(WA[2,1]*x + WA[2,2]));
return HyperellipticCurve(g / d);

end intrinsic;


intrinsic FactorReconstructG2(P::., K::Fld) -> .
{Reconstructs genus 2 factor.}

// Recover small period matrix:
CC := Parent(P[1,1]); RCC<xCC> := PolynomialRing(CC);
g := #Rows(Transpose(P));
Omega1 := Submatrix(P, 1, 1, g, g); Omega2 := Submatrix(P, g + 1, 1, g, g);
PSmall := Omega2 * Omega1^(-1);

rosensCC := RosenhainInvariants(PSmall);
fCC := xCC * (xCC - 1) * &*[ xCC - rosenCC : rosenCC in rosensCC ];
g2sCC := G2Invariants(HyperellipticCurve(fCC));
g2s := [ AlgebraizeElementInRelativeField(g2CC, K) : g2CC in g2sCC ];
X := HyperellipticCurveFromG2Invariants(g2s);
X := TwistDifferentialBasis(X, P);
return X;

end intrinsic;


intrinsic FactorsFromProjections(projs::List) -> List
{Recovers algebraic expressions for the factors.}

Facs := [* *];
for proj in projs do
    lat := proj[1]; K := BaseRing(proj[2][1]);
    g := #Rows(Transpose(lat));
    if g eq 1 then
        Append(~Facs, FactorReconstructG1(lat, K));
    elif g eq 2 then
        Append(~Facs, FactorReconstructG2(lat, K));
    else
        Append(~Facs, "");
    end if;
end for;
return Facs;

end intrinsic;
