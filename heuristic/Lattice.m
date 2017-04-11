/**
 *  Lattices of subfield data
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
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


intrinsic EndomorphismLattice(GeoEndList::List) -> List
{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}

AsAlg, Rs, As := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);
Gp, Gf, Gphi := AutomorphismGroup(L);

Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~Hs, CompareGroups); Reverse(~Hs);

LatReps := [* *]; LatStructs := [* *]; LatDescs := [* *];
for H in Hs do
    Hf := [ Gphi(genH) : genH in Generators(H) ];
    LatRep := [* *]; LatStruct := [* *]; LatDesc := [* *];
    K := FixedField(L, Hf);
    // TODO: Indicate class group and treat the relative case (scaffolding in place).
    K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
    Append(~LatRep, K); Append(~LatStruct, K); Append(~LatDesc, K_seq);
    EndoReps := EndomorphismBasis(GeoEndList, Hf);
    EndoStructs, EndoDescs := EndomorphismStructure(EndoReps);
    LatRep cat:= EndoReps; LatStruct cat:= EndoStructs; LatDesc cat:= EndoDescs;
    SatoTate := SatoTateGroup(GeoEndList, Hf);
    Append(~LatRep, SatoTate); Append(~LatStruct, SatoTate); Append(~LatDesc, SatoTate);
    Append(~LatReps, LatRep); Append(~LatStructs, LatStruct); Append(~LatDescs, LatDesc);
end for;

return [* LatReps, LatStructs, LatDescs *];

end intrinsic;
