/***
 *  Lattice functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function EndomorphismLatticeG2SingleElement(K, AsAlg, As, Rs : AddTensor :=
    true, AddRing := true, AddSatoTate := true, AddDecomposition := true);
// Input:   Two sets of matrices representing a geometric endomorphism basis
//          analytically and rationally.
// Output:  A single element of the full lattice of endomorphism data, labelled
//          by the corresponding subgroup and subfield, along with non-trivial
//          idempotents.

GeoEndList := [* AsAlg, As, Rs *];
L := BaseRing(AsAlg[1]);
EndRep, idemsList, GeoFactorsQQ := EndomorphismData(GeoEndList, L, L, 0, 0, 0,
    0, 0 : AddTensor := true, AddRing := AddRing, AddSatoTate := AddSatoTate,
    AddDecomposition := AddDecomposition);
Shorthand := GeoEndRRShorthand(EndRep[2]);
if Degree(K) eq 1 then
    Gp, Gf, Gphi := AutomorphismGroup(L);
    GensHp := Generators(Gp);
    KL := sub<L | 1>;
    return EndomorphismData(GeoEndList, L, KL, Gp, GensHp, Gphi, GeoFactorsQQ,
        Shorthand : AddTensor := AddTensor, AddRing := AddRing, AddSatoTate :=
        AddSatoTate, AddDecomposition := AddDecomposition);
end if;
// Take smallest subfield of K that fits inside L and take corresponding fixed
// group
S := Subfields(K);
S := [ tup[1] : tup in S ];
Sort(~S, CompareSubfield);
Reverse(~S);
for M in S do
    test, f := IsSubfield(M, L);
    if test then
        KL := sub<L | f(M.1)>;
        break;
    end if;
end for;
Gp, Gf, Gphi := AutomorphismGroup(L);
if test then
    //GensHp := Generators(FixedGroupFixed(L, KL, Gp, Gphi));
    GensHp := Generators(FixedGroup(L, KL));
else
    KL := sub<L | 1>;
    GensHp := Generators(Gp);
end if;

return EndomorphismData(GeoEndList, L, KL, Gp, GensHp, Gphi, GeoFactorsQQ,
    Shorthand : AddTensor := AddTensor, AddRing := AddRing, AddSatoTate :=
    AddSatoTate, AddDecomposition := AddDecomposition);

end function;


function EndomorphismLatticeG2(AsAlg, As, Rs : Geometric := true, AddTensor :=
    false, AddRing := false, AddSatoTate := false, AddDecomposition := false);
// Input:   Two sets of matrices representing a geometric endomorphism basis
//          analytically and rationally.
// Output:  The full lattice of endomorphism data, labelled by the corresponding
//          subgroup and subfield, along with non-trivial idempotents.
//          Fields of definition for the latter is also returned.

GeoEndList := [* AsAlg, As, Rs *];
L := BaseRing(AsAlg[1]);
// Simplest case (we get around the lacking string because it gets calculated
// during the run) over the field of definition of the endomorphism ring:
// TODO: Note that if we do not run this with AddTensor set, we do not get
// uniform output. I do not want to bother with this right now.
EndRep, idemsList, GeoFactorsQQ := EndomorphismData(GeoEndList, L, L, 0, 0, 0,
    0, 0 : AddTensor := true, AddRing := AddRing, AddSatoTate := AddSatoTate,
    AddDecomposition := AddDecomposition);
if Geometric or (Degree(L) eq 1) then
    return [* [* Eltseq(L.1) *] cat EndRep *], Eltseq(L.1), idemsList;
end if;
Shorthand := GeoEndRRShorthand(EndRep[2]);
// Two heavy operations follow:
Gp, Gf, Gphi := AutomorphismGroup(L);
Ktups := Seqlist(Subfields(L));
// Have to add in the rationals manually
// FIXME: Magma should really add it...
QQ := RationalsAsNumberField();
iotaQQ := EmbeddingMap(QQ, L);
Append(~Ktups, <QQ, iotaQQ>);
EndReps := [* [* Eltseq(L.1) *] cat EndRep *];
idemsList0 := idemsList;
K0 := L;
// We construct a list of lists EndReps. Its entries should be thought of as the
// subfields of the field of definition L of the full endomorphism ring,
// along with all the endomorphisms defined over those subfields.
for tup in Ktups do
    K := tup[1];
    if K ne L then
        // TODO: Check Magma fix
        //Hp := FixedGroupFixed(L, K, Gp, Gphi);
        Hp := FixedGroup(L, K);
        GensHp := Generators(Hp);
        // TODO: Add description of output of upcoming function
        EndRep, idemsList := EndomorphismData(GeoEndList, L, K, Gp, GensHp,
            Gphi, GeoFactorsQQ, Shorthand : AddTensor := AddTensor, AddRing :=
            AddRing, AddSatoTate := AddSatoTate, AddDecomposition := AddDecomposition);
        // Keep track of field of definition of decomposition so as to return it
        // over the smallest possible field later:
        if #idemsList[1] ne 0 and Degree(K) le Degree(K0) then
            idemsList0 := idemsList;
            K0 := K;
        end if;
        Append(~EndReps, [* Eltseq(L ! K.1) *] cat EndRep);
    end if;
end for;

// We later convert this, but for now we return "structured" data:
return EndReps, Eltseq(L ! K0.1), idemsList0;

end function;

