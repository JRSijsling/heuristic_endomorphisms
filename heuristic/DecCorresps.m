/***
 *  Correspondences from factors
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


// TODO: Make singular, use new data structure indices
// Export first

intrinsic CorrespondencesFromFactorsAndProjections(X::Crv, facs::List, projs::List) -> List
{Determines morphisms from factors and analytic projections.}

Mors := [* *];
for i:=1 to #facs do
    fac := facs[i]; proj := projs[i];
    lat := proj[1]; g := #Rows(Transpose(lat));
    if g eq 1 then
        Append(~Mors, CorrespondenceG1(X, fac, proj));
    elif g eq 2 then
        Append(~Mors, CorrespondenceGeneral(X, fac, proj));
    else
        Append(~Mors, "");
    end if;
end for;
return true, Mors;

end intrinsic;


intrinsic CorrespondenceG1(X::Crv, Y::Crv, proj::List) -> .
{Determines morphism and verifies in genus 1.}

A := proj[2][1] * IdempotentDenominator(proj[2][2]);
XL, AsL, PL := NonWeierstrassBasePoint(X, [ A ]); AL := AsL[1];
L := BaseRing(AL); YL := ChangeRing(Y, L); QL := YL ! [1, 0, 0];
return CantorMorphismFromMatrixSplit(XL, PL, YL, QL, AL);

end intrinsic;


intrinsic CorrespondenceGeneral(X::Crv, Y::Crv, proj::List) -> .
{Determines morphism and verifies in general.}

// TODO: FindPoint needs improvement

return "";

end intrinsic;
