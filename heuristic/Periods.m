/***
 *  Determining period matrices
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

/* Enable Oldenburg if you have access to the relevant code by Pascal Molin,
 * Christian Neurohr et al. */

intrinsic PeriodMatrixHyperelliptic(f::RngUPolElt, h::RngUPolElt : HaveOldenburg := false) -> AlgMatElt
{Computes a (big) period matrix of the hyperelliptic curve defined by f and h.}

g := 4*f + h^2;
if not HaveOldenburg then
    J := AnalyticJacobian(g);
    return Transpose(BigPeriodMatrix(J));
end if;
return Transpose(PeriodMatrix(g));

end intrinsic;


intrinsic PeriodMatrixPlane(F::RngMPolElt : HaveOldenburg := false) -> AlgMatElt
{Computes a (big) period matrix of the plane curve defined by F.}

SCC<x0,x1,x2> := Parent(F); CC := BaseRing(SCC); RCC<x,y> := PolynomialRing(CC, 2);
h := hom<SCC -> RCC | [x,y,1]>; f := h(F);
if not HaveOldenburg then
    return 0;
end if;
return Transpose(PeriodMatrix(f));

end intrinsic;
