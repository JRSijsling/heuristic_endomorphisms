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


intrinsic EndomorphismLattice(GeoEndoRep::SeqEnum) -> SeqEnum
{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}

gensTan := [ gen[1] : gen in GeoEndoRep ];
gensHom := [ gen[2] : gen in GeoEndoRep ];
gensApp := [ gen[3] : gen in GeoEndoRep ];
L := BaseRing(gensTan[1]);
Gp, Gf, Gphi := AutomorphismGroup(L);

// TODO: Reverse or not (we want to calculate the geometric case first and can always reverse later)
Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~Hs, CompareGroups); Reverse(~Hs);

Lat := [ ];
for H in Hs do
    OverK := [* *];
    gensH := Generators(H); GalK := [* gensH, Gphi *];
    K := FixedField(L, [ Gphi(genH) : genH in gensH ]);
    // TODO: Indicate class group and treat the relative case (scaffolding in place).
    K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
    K_desc := [* K_seq, K *];
    Append(~OverK, K_desc);
    Append(~OverK, EndomorphismStructure(GeoEndoRep, GalK));
    Append(~Lat, OverK);
end for;
return Lat;

end intrinsic;
