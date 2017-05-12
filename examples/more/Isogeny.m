load "Initialize.m";
load "heuristic/Isogeny.m";

prec0 := 100;
CC<I> := ComplexField(prec0);
R<x> := PolynomialRing(CC);
fC := x^4 + x^2;
hC := x^3 + 1;
gC := fC + (hC/2)^2;
C := HyperellipticCurve(gC);
JC := AnalyticJacobian(C);
PC := BigPeriodMatrix(JC);
PC := Transpose(PC);
//print "Period matrix of C:";
//print PC;

fE := x^3 - x;
hE := x + 1;
gE := fE + (hE/2)^2;
E := HyperellipticCurve(gE);
JE := AnalyticJacobian(E);
PE := BigPeriodMatrix(JE);
PE := Transpose(PE);
//print "Period matrix of E:";
//print PE;

fE := -2*x^4 + 4*x^2 - 9*x - 14;
hE := x^3 + 1;
gE := fE + (hE/2)^2;
E := HyperellipticCurve(gE);
JE := AnalyticJacobian(E);
PE := BigPeriodMatrix(JE);
PE := Transpose(PE);
//print "Period matrix of E:";
//print PE;

JC := ComplexStructure(PC);
JE := ComplexStructure(PE);
M := RationalEndomorphismEquations(JC, JE);

/*
RR := RealField(CC);
M := Matrix(RR, [ [ Sqrt(RR ! 2) ], [ 2 * Sqrt(RR ! 2) ] ]);
print "Kernel:";
print IntegralKernel(M);
*/

K := IntegralKernel(M);
print K;

RR := RealField(CC);
v := K[1];
//M := Matrix(RR, 2, 4, Eltseq(v));
//M := Matrix(RR, 4, 2, Eltseq(v));
M := Matrix(RR, 4, 4, Eltseq(v));
print "Small?";
//print M*JC - JE*M;
print M*JE - JC*M;

exit;
