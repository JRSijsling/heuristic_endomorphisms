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

CC := BaseRing(Parent(F));
return Transpose(PeriodMatrix(F));

end function;
