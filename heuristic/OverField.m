/***
 *  Subrings of the geometric endomorphism ring over subfields
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


intrinsic EndomorphismBasis(GeoEndList::List, GensHf::SeqEnum) -> List
{Extracts basis over subfield determines by a list of automorphisms.}

// Setting the stage
AsAlg, As, Rs := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);
N := #AsAlg;

// Trivial case
if #GensHf eq 0 then
    return GeoEndList;
end if;

// The vector space representing the full endomorphism algebra
Kern := VectorSpace(Rationals(), N);
// Successively filter by writing down the conditions for a matrix to be fixed
// under a given generator
for sigma in GensHf do
    Msigma := Matrix([MatrixInBasis(ConjugateMatrix(sigma, A), AsAlg) : A in AsAlg]);
    Msigma := Msigma - IdentityMatrix(Rationals(), N);
    Kern := Kern meet Kernel(Msigma);
end for;

// Retrieve representations for a basis of of the endomorphism ring by taking a
// pure lattice / saturation
Lat := PureLattice(Lattice(Matrix(Basis(Kern))));
B := Basis(Lat);

// Constructing said basis
AsAlg0 := [&+[b[n] * AsAlg[n] : n in [1..N]] : b in B];
As0 := [&+[b[n] * As[n] : n in [1..N]] : b in B];
Rs0 := [&+[b[n] * Rs[n] : n in [1..N]] : b in B];
return [* AsAlg0, As0, Rs0 *];

end intrinsic;


function CompareFields(K1, K2);
// Input:   Two subfields, fields, or polynomials.
// Output:  A comparison function: field with smaller degrees are smaller.

if Degree(K1) lt Degree(K2) then
    return -1;
elif Degree(K1) eq Degree(K2) then
    return 0;
else
    return 1;
end if;

end function;


intrinsic EndomorphismBasis(GeoEndList::List, K::Fld) -> List
{Extracts basis over a general field.}

// Setting the stage
AsAlg, As, Rs := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);

// Trivial case
if Degree(K) eq 1 then
    Gp, Gf, Gphi := AutomorphismGroup(L);
    GensHf := [ Gphi(gen) : gen in Generators(Gp) ];
    return EndomorphismBasisOverSubfield(GensHf, GeoEndList);
end if;

// Take smallest subfield of K that fits inside L and take corresponding fixed group
S := Subfields(K); S := [ tup[1] : tup in S ];
Sort(~S, CompareSubfield); Reverse(~S);
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
    // Case where we need to go down all the way
    GensHf := [ Gphi(gen) : gen in Generators(Gp) ];
end if;
return EndomorphismBasis(GensHf, GeoEndList);

end intrinsic;
