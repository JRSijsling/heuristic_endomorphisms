/***
 *  Algebraizing the basis over a field
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function EndomorphismBasisOverSubfield(GensHf, GeoEndList);
// Input:   A set of generators of a Galois group and a sequence of matrices
//          that represent the geometric endomorphism basis. The Galois group in
//          question is a subgroup of that of the field of definition of the
//          entries of AsAlg.
// Output:  A basis of the endomorphism ring that is defined over the fixed
//          field of the given generators, again a sequence.

// Trivial case first:
if #GensHf eq 0 then
    return GeoEndList;
end if;
AsAlg, As, Rs := Explode(GeoEndList);
N := #AsAlg;
L := BaseRing(AsAlg[1]);
// The vector space representing the full endomorphism algebra:
Kern := VectorSpace(Rationals(), N);
// Successively filter by writing down the conditions for a matrix to be fixed
// under a given generator:
for sigma in GensHf do
    Msigma := Matrix([MatrixRatInBasisOverNF(ConjugateMatrix(sigma, A), AsAlg) :
        A in AsAlg]);
    Msigma := Msigma - IdentityMatrix(Rationals(), N);
    Kern := Kern meet Kernel(Msigma);
end for;
// Retrieve representations for a basis of of the endomorphism ring by taking a
// pure lattice / saturation:
Lat := PureLattice(Lattice(Matrix(Basis(Kern))));
B := Basis(Lat);
// Constructing said basis:
AsAlg0 := [&+[b[n] * AsAlg[n] : n in [1..N]] : b in B];
As0 := [&+[b[n] * As[n] : n in [1..N]] : b in B];
Rs0 := [&+[b[n] * Rs[n] : n in [1..N]] : b in B];

return [* AsAlg0, As0, Rs0 *];

end function;


function EndomorphismBasisOverField(K, AsAlg, As, Rs);
// Input:   A field and a field of definition of a geometric endomorphism basis.
// Output:  A basis of the endomorphism ring that is defined over the fixed
//          field of the given field, again a sequence.

L := BaseRing(AsAlg[1]);
GeoEndList := [* AsAlg, As, Rs *];
if Degree(K) eq 1 then
    Gp, Gf, Gphi := AutomorphismGroup(L);
    GensHf := [ Gphi(gen) : gen in Generators(Gp) ];
    return EndomorphismBasisOverSubfield(GensHf, GeoEndList);
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
    GensHf := [ Gphi(gen) : gen in Generators(FixedGroup(L, KL)) ];
else
    GensHf := [ Gphi(gen) : gen in Generators(Gp) ];
end if;

return EndomorphismBasisOverSubfield(GensHf, GeoEndList);

end function;
