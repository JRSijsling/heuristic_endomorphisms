t := 2;
t := 2/3;
R<x> := PolynomialRing(Rationals());
F<s> := NumberField(x^2 + 3 - 14*t^2 + 27*t^4);
R<x> := PolynomialRing(F);

A := s / (2*t);
B := (1 + 3*t^2)/(1 - 3*t^2);
Q := -(1 - 2*t^2 + 9*t^4)*(1 - 28*t^2 + 166*t^4 - 252*t^6 + 81*t^8)/(4*t^2*(1 - 3*t^2)^2*(1 - t^2)*(1 - 9*t^2));
f := x*(x^4 + (A - B)*x^3 + Q*x^2 + (A + B)*x + 1);

X := HyperellipticCurve(f);
I := G2Invariants(X);
I := [ Rationals() ! inv : inv in I ];
print I;
//print HyperellipticCurveFromG2Invariants(I);

R<x> := PolynomialRing(Rationals());
F<r> := NumberField(x^2 + 1);
A<x,y> := AffineSpace(F, 2);
f := y^2 + 3 - 14*x^2 + 27*x^4;
C := Curve(Scheme(A, f));
C := ProjectiveClosure(C);
P := C ! [1, -4*r, 1];
E, phi := EllipticCurve(C, P);

U,f := TorsionSubgroup(E);
u_exc := [ 0*U.1, U.2 ];
for u in U do
    print u;
    if not u in u_exc then
        Ps := Points(f(u) @@ phi);
        for P in Ps do
            R<x> := PolynomialRing(F);
            s := P[1]; t := P[2];
            if 4*t^2*(1 - 3*t^2)^2*(1 - t^2)*(1 - 9*t^2) ne 0 then
                print P;

                A := s / (2*t);
                B := (1 + 3*t^2)/(1 - 3*t^2);
                Q := -(1 - 2*t^2 + 9*t^4)*(1 - 28*t^2 + 166*t^4 - 252*t^6 + 81*t^8)/(4*t^2*(1 - 3*t^2)^2*(1 - t^2)*(1 - 9*t^2));
                p := x*(x^4 + (A - B)*x^3 + Q*x^2 + (A + B)*x + 1);

                X := HyperellipticCurve(p);
                I := G2Invariants(X);
                I := [ Rationals() ! inv : inv in I ];
                print I;
            end if;
        end for;
    end if;
end for;

