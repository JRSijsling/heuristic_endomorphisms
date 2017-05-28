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


intrinsic EndomorphismLattice(GeoEndoRep::SeqEnum : Optimize := false) -> List
{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}

gensTan := [ gen[1] : gen in GeoEndoRep ];
gensHom := [ gen[2] : gen in GeoEndoRep ];
gensApp := [ gen[3] : gen in GeoEndoRep ];
L<s> := BaseRing(gensTan[1]); F := BaseField(L);
Gp, Gf, Gphi := AutomorphismGroup(L);

Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~Hs, CompareGroups);

entries := [ ];
// The code of this first (geometric) step is a copy of that below except for
// shorthand extraction. Of course it could be simpler, but there is no time
// loss as a result.
entry := [* *];
gensH := [ ]; GalK := [* gensH, Gphi *];
K<s> := L;

if IsRational(F) then
    K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
else
    K_seq := [ [ Integers() ! c : c in Eltseq(cs) ] : cs in Eltseq(MinimalPolynomial(K.1)) ];
end if;
K_desc := [* K_seq, K *];
EndoStruct := EndomorphismStructure(GeoEndoRep, GalK);
Append(~entry, K_desc);
Append(~entry, EndoStruct);
//Append(~entry, ClassNumber(AbsoluteField(K)));
Append(~entries, entry);
Shorthand := SatoTateShorthand(EndoStruct);

for H in Hs[2..#Hs] do
    entry := [* *];
    gensH := Generators(H); GalK := [* gensH, Gphi *];
    if HasRationalBase(L) then
        K := MakeRelative(FixedField(L, [ Gphi(genH) : genH in gensH ]), Rationals());
    else
        K := RelativeFixedField(L, [ Gphi(genH) : genH in gensH ]);
    end if;

    K := ClearFieldDenominator(K);
    if (not IsRational(K)) and Optimize then
        K := OptimizedRepresentation(K);
        K := ClearFieldDenominator(K);
    end if;
    K<s> := K;

    if IsRational(F) then
        K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
    else
        K_seq := [ [ Integers() ! c : c in Eltseq(F ! coeff) ] : coeff in Eltseq(MinimalPolynomial(K.1)) ];
    end if;

    K_desc := [* K_seq, K *];
    EndoStruct := EndomorphismStructure(GeoEndoRep, GalK : Shorthand := Shorthand);
    Append(~entry, K_desc);
    Append(~entry, EndoStruct);
    //Append(~entry, ClassNumber(AbsoluteField(K)));
    Append(~entries, entry);
end for;

F_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(F.1)) ];
base := [* F_seq, F *];
return [* base, entries *];

end intrinsic;
