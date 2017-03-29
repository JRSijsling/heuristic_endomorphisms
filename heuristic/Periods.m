function PeriodMatrixHyperelliptic(f, h)
// Input:   Two polynomials f and h over the rationals, reals or CC.
// Output:  The corresponding period matrix.

g := 4*f + h^2;
// NOTE: Comment out next line if you do not have the Oldenburg functionality
return Transpose(PeriodMatrix(g));
J := AnalyticJacobian(g);
return Transpose(BigPeriodMatrix(J));

end function;


function PeriodMatrixPlane(F)
// Input:   A polynomials over the rationals, reals or CC.
// Output:  The corresponding period matrix.

SCC<x0,x1,x2> := Parent(F);
CC := BaseRing(SCC);
RCC<x,y> := PolynomialRing(CC, 2);
h := hom<SCC -> RCC | [x,y,1]>;
f := h(F);
return Transpose(PeriodMatrix(f));

end function;
