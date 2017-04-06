/***
 *  Examples of genus 2 factors
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

AttachSpec("../spec");

prec := 300;
CC := ComplexFieldExtra(prec);
RCC := PolynomialRing(CC);
R<x> := PolynomialRing(QQ);
f := x^5 + x^4 - x^3 + 3*x^2 + 4*x + 1;
f := 2*x^6 - 6*x^5 + x^4 - 3*x^3 + 3*x^2 + 4*x + 7;
fCC := RCC ! f;

XCC := HyperellipticCurve(fCC);
JCC := AnalyticJacobian(XCC);
PBig := BigPeriodMatrix(JCC);
PSmall := SmallPeriodMatrix(JCC);

result := HyperellipticCurveFromBigPeriodMatrixG2(PBig, PSmall);
print f;
print result;

exit;
