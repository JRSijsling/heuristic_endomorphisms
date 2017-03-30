function RecognizeG2Invariants(g2sCC : epscomp := epscomp0, epsLLL := epscomp0);
// Recognizes a bunch of G2Invariants, together, as algebraic elements

pols := [ PolynomializeElement(g2 : epscomp := epscomp, epsLLL := epsLLL) : g2 in g2sCC ];
K := MySplittingField(pols);
f := MinimalPolynomial(K.1); frep := Eltseq(f);
CC := Parent(g2sCC[1]); h := Roots(f, CC)[1][1];
g2s := [ AlgebraizeElementInField(g2CC, frep, h : epscomp := epscomp, epsLLL := epsLLL) : g2CC in g2sCC ];
if Degree(K) eq 1 then
    g2s := [ Rationals() ! g2 : g2 in g2s ];
end if;
return g2s;

end function;


function TwistingFactor(P, X);
// TODO: Just generic for now

CC := Parent(P[1,1]); RCC := PolynomialRing(CC);
f, h := HyperellipticPolynomials(X);
Q := Transpose(BigPeriodMatrix(AnalyticJacobian(HyperellipticCurve(RCC ! f))));
Q := Matrix(CC, #Rows(Q), #Rows(Transpose(Q)), Eltseq(Q));
As, Rs := GeometricIsogenyBasisFromPeriodMatrices(P, Q);
Apol := PolynomializeMatrix(As[1]);
K := MySplittingField(Eltseq(Apol));
return SquareFree(Discriminant(Integers(K)));

end function;


function HyperellipticCurveFromBigPeriodMatrixG2(PBig, PSmall : epscomp := epscomp0, epsLLL := epscomp0);
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
g2s := RecognizeG2Invariants(g2sCC : epscomp := epscomp, epsLLL := epsLLL);
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

end function;

