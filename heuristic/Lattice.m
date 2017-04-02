/**
 *  Lattice functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function EndomorphismLatticeG3SingleElement(K, AsAlg, As, Rs)
// Input:   Two sets of matrices representing a geometric endomorphism basis
//          analytically and rationally.
// Output:  A single element of the full lattice of endomorphism data, labelled
//          by the corresponding subgroup and subfield, along with non-trivial
//          idempotents.

// TODO: Make relative

GeoEndList := [* AsAlg, As, Rs *];
L := BaseRing(AsAlg[1]);
if Degree(L) eq 1 then
    return EndomorphismDataG3(GeoEndList, L, L);
end if;
// Take smallest subfield of K that fits inside L and take corresponding fixed
// group
S := Subfields(K); S := [ tup[1] : tup in S ];
Sort(~S, CompareSubfield); Reverse(~S);
for M in S do
    test, f := IsSubfield(M, L);
    if test then
        KL := sub<L | f(M.1)>;
        break;
    end if;
end for;
if not test then
    // Trivial subextension:
    KL := sub<L | 1>;
end if;

return EndomorphismDataG3(GeoEndList, L, KL);

end function;


function EndomorphismLatticeG3(AsAlg, As, Rs);
// Input:   Two sets of matrices representing a geometric endomorphism basis
//          analytically and rationally.
// Output:  The full lattice of endomorphism data, labelled by the corresponding
//          subgroup and subfield, along with non-trivial idempotents.
//          Fields of definition for the latter is also returned.

GeoEndList := [* AsAlg, As, Rs *];
L := BaseRing(AsAlg[1]);
// Simplest case (we get around the lacking string because it gets calculated
// during the run) over the field of definition of the endomorphism ring:
if Degree(L) eq 1 then
    EndRep := EndomorphismDataG3(GeoEndList, L, L);
    return [* [* Eltseq(L.1) *] cat EndRep *];
end if;
// Two heavy operations follow:
Gp, Gf, Gphi := AutomorphismGroup(L);
Ktups := Seqlist(Subfields(L));
// Have to add in the rationals manually
// FIXME: Magma should really add it...
QQ := RationalsAsNumberField(); iotaQQ := EmbeddingMap(QQ, L);
Append(~Ktups, <QQ, iotaQQ>);
// We construct a list of lists EndReps. Its entries should be thought of as the
// subfields of the field of definition L of the full endomorphism ring,
// along with all the endomorphisms defined over those subfields.
EndRep := EndomorphismDataG3(GeoEndList, L, L);
EndReps := [* [* Eltseq(L.1) *] cat EndRep *];
for tup in Ktups do
    K := tup[1];
    if K ne L then
        EndRep := EndomorphismDataG3(GeoEndList, L, K);
        Append(~EndReps, [* Eltseq(L ! K.1) *] cat EndRep);
    end if;
end for;

// We later convert this, but for now we return "structured" data:
return EndReps;

end function;
