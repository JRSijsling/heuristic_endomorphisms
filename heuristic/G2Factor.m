/***
 *  Genus 2 curve from decomposition
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function RecognizeG2Invariants(g2sCC)
// Recognizes a bunch of G2Invariants, together, as algebraic elements

CC := Parent(g2sCC[1]);
QQ := Rationals(); SetInfinitePlace(QQ, InfinitePlaces(QQ)[1]);
pols := [ RelativeMinimalPolynomial(g2, QQ) : g2 in g2sCC ];
K := RelativeSplittingFieldExtra(pols);
g2s := [ AlgebraizeElementInRelativeField(g2CC, K) : g2CC in g2sCC ];
return g2s;

end function;


function TwistingFactor(P, X)
// TODO: Just generic for now

CC := Parent(P[1,1]); RCC := PolynomialRing(CC);
QQ := Rationals(); SetInfinitePlace(QQ, InfinitePlaces(QQ)[1]);
f, h := HyperellipticPolynomials(X);
Q := Transpose(BigPeriodMatrix(AnalyticJacobian(HyperellipticCurve(RCC ! f))));
Q := Matrix(CC, #Rows(Q), #Rows(Transpose(Q)), Eltseq(Q));
AsRs := GeometricIsogenyBasisApproximations(P, Q);
pols := RelativeMinimalPolynomialsMatrices(AsRs[1], QQ);
K := SplittingField(pols);
return SquareFree(Discriminant(Integers(K)));

end function;


intrinsic HyperellipticCurveFromBigPeriodMatrixG2(PBig::., PSmall::.) -> .
{Finds a genus 2 factor.}
// (P::AlgMatElt: prec := 300) -> CrvHyp

PBig := Transpose(PBig);
//g := #Rows(Transpose(PBig));
//Omega1 := Submatrix(PBig, 1, 1, g, g);
//Omega2 := Submatrix(PBig, g + 1, 1, g, g);
//PSmall := Omega2 * Omega1^(-1);
rosensCC := RosenhainInvariants(PSmall);
CC := Parent(PSmall[1,1]);
RCC<xCC> := PolynomialRing(CC);
fCC := xCC * (xCC - 1) * &*[ xCC - rosenCC : rosenCC in rosensCC ];
g2sCC := G2Invariants(HyperellipticCurve(fCC));
g2s := RecognizeG2Invariants(g2sCC);
X := HyperellipticCurveFromG2Invariants(g2s);
if BaseRing(X) eq Rationals() then
    X := ReducedMinimalWeierstrassModel(X);
    f, h := HyperellipticPolynomials(X);
    f /:= LeadingCoefficient(f);
    X := HyperellipticCurve(f);
    d := TwistingFactor(PBig, X);
    X := HyperellipticCurve(d*f);
end if;
return X;

end intrinsic;
