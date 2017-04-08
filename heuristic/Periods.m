/***
 *  Determining period matrices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


/* Enable Oldenburg if you have access to the relevant code by Pascal Molin,
 * Christian Neurohr et al. */

intrinsic PeriodMatrix(f::RngUPolElt, h::RngUPolElt : HaveOldenburg := false) -> AlgMatElt
{Computes a (big) period matrix of the hyperelliptic curve defined by f and h.}

g := 4*f + h^2;
if not HaveOldenburg then
    J := AnalyticJacobian(g);
    return Transpose(Matrix(BaseRing(g), BigPeriodMatrix(J)));
end if;
return Transpose(Matrix(BaseRing(g), PeriodMatrix(g)));

end intrinsic;


intrinsic PeriodMatrix(F::RngMPolElt : HaveOldenburg := false) -> AlgMatElt
{Computes a (big) period matrix of the plane curve defined by F.}

SCC<x0,x1,x2> := Parent(F); CC := BaseRing(SCC); RCC<x,y> := PolynomialRing(CC, 2);
h := hom<SCC -> RCC | [x,y,1]>; f := h(F);
if not HaveOldenburg then
    return 0;
end if;
return Transpose(Matrix(BaseRing(F), PeriodMatrix(f)));

end intrinsic;
