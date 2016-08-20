R<t> := PolynomialRing(Rationals());
K<r> := NumberField(t^2 + 55);

R2<x,y> := PolynomialRing(K, 2);
f1 := (27/2540*r*x + 81/2540*r)*y + (-15/508*x^4 + 279/2540*x^3 + 963/5080*x^2 - 2943/2540*x - 891/508);
f2 := (-24129/161290*x^5 + 106893/1290320*x^4 + 46672821/12903200*x^3 + 60865587/12903200*x^2 - 53804007/6451600*x - 10381203/806450)*y + (-567/32258*r*x^8 - 63477/6451600*r*x^7 - 802521/12903200*r*x^6 - 9287703/12903200*r*x^5 + 7077051/6451600*r*x^4 + 112291029/12903200*r*x^3 + 80322273/12903200*r*x^2 - 78102873/6451600*r*x - 10381203/806450*r);
f3 := x^4 + 171/127*x^3 - 26871/5080*x^2 - 15201/1270*x - 8901/1270;

PC<XC,YC,ZC> := ProjectiveSpace(K, [1,3,1]);
RC<XC,YC,ZC> := CoordinateRing(PC);

PE<XE,YE,ZE> := ProjectiveSpace(K, [1,2,1]);
RE<XE,YE,ZE> := CoordinateRing(PE);

PP<XP,ZP> := ProjectiveSpace(K, [1,1]);
RP<XP,ZP> := CoordinateRing(PP);

pC := x^6 - 8*x^4 + 2*x^3 + 16*x^2 - 36*x - 55;
pE := -10582/27*x^4 + 215/3*x^3 + x;

FC := YC^2 - RC ! (ZC^6*Evaluate(pC, [XC/ZC, YC/ZC^2]));
FE := YE^2 - RE ! (ZE^4*Evaluate(pE, [XE/ZE, YE/ZE^2]));

C := Curve(PC, FC);
E := Curve(PE, FE);
P := Curve(PP);

F1 := RC ! (ZC^4*Evaluate(f1, [XC/ZC, YC/ZC^3]));
F2 := RC ! (ZC^8*Evaluate(f2, [XC/ZC, YC/ZC^3]));
F3 := RC ! (ZC^4*Evaluate(f3, [XC/ZC, YC/ZC^3]));
print F1;
print F2;
print F3;

h1 := map<C -> E | [F1, F2, F3]>;
print Degree(h1);

h2 := map<C -> P | [F1, F3]>;
print Degree(h2);
print h2;

exit;
