/**
 *  Lattices of subfield data
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


function CompareGroups(G1, G2);
// Input:   Two subgroups or groups.
// Output:  A comparison function: groups with smaller cardinality are smaller.

if #G1 lt #G2 then
    return -1;
elif #G1 eq #G2 then
    return 0;
else
    return 1;
end if;

end function;


intrinsic EndomorphismLatticeDescription(GeoEndList::List) -> List
{Returns the lattice of descriptions per subfield.}

AsAlg, As, Rs := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);
Gp, Gf, Gphi := AutomorphismGroup(L);
Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~S, CompareGroups); Reverse(~S);

Lat := [* *];
for H in Hs do
    K := FixedField(H);
    EndoDescList := EndomorphismData(EndomorphismBasis(GeoEndList, K)[2]);
    Append(~Lat, [* Eltseq(L ! K.1) *] cat EndoDescList);
end for;

return Lat;

end intrinsic;
