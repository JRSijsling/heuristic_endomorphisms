"""
 *  Some examples of endomorphism rings of genus 2 hyperelliptic curves
 *
 *  Copyright (C) 2016, 2017 Edgar Costa,
 *                           Davide Lombardo,
 *                           Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
"""

# Add if no initilization script set:
load('Initialize.sage')

prec = 300
CC = magma.ComplexField(prec)
RCC = magma.PolynomialRing(CC)
R.<x> = PolynomialRing(QQ)
f = x^5 + x^4 - x^3 + 3*x^2 + 4*x + 1
f = 2*x^6 - 6*x^5 + x^4 - 3*x^3 + 3*x^2 + 4*x + 7
fCC = RCC(f)
XCC = magma.HyperellipticCurve(fCC)
JCC = magma.AnalyticJacobian(XCC)
PBig = magma.BigPeriodMatrix(JCC)
PSmall = magma.SmallPeriodMatrix(JCC)
result = magma.HyperellipticCurveFromBigPeriodMatrixG2(PBig, PSmall);
print result

