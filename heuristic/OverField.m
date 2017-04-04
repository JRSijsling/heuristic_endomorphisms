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
AsAlg, Rs, As := Explode(GeoEndList);
L := BaseRing(AsAlg[1]); K := FixedField(L, GensHf);
N := #AsAlg;

// Trivial case
if #GensHf eq 0 then
    return GeoEndList;
end if;

// The vector space representing the full endomorphism algebra
Ker := VectorSpace(Rationals(), N);
// Successively filter by writing down the conditions for a matrix to be fixed
// under a given generator
for sigma in GensHf do
    Msigma := Matrix([MatrixInBasis(ConjugateMatrix(sigma, A), AsAlg) : A in AsAlg]);
    Msigma -:= IdentityMatrix(Rationals(), N);
    Ker meet:= Kernel(Msigma);
end for;

// Retrieve representations for a basis of of the endomorphism ring by taking a
// pure lattice / saturation
Lat := PureLattice(Lattice(Matrix(Basis(Ker))));
B := Basis(Lat);

// Constructing said basis
AsAlg := [ &+[ b[n] * AsAlg[n] : n in [1..N] ] : b in B ];
Rs    := [ &+[ b[n] * Rs[n]    : n in [1..N] ] : b in B ];
As    := [ &+[ b[n] * As[n]    : n in [1..N] ] : b in B ];
AsAlg := [ Matrix(K, A) : A in AsAlg ];
return [* AsAlg, Rs, As *];

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
if (Degree(K) eq 1) or (Type(L) eq FldRat) then
    Gp, Gf, Gphi := AutomorphismGroup(L);
    GensHf := [ Gphi(gen) : gen in Generators(Gp) ];
    return EndomorphismBasis(GeoEndList, GensHf);
end if;

// Take smallest subfield of K that fits inside L and take corresponding fixed group
S := Subfields(K); S := [ tup[1] : tup in S ];
Sort(~S, CompareFields); Reverse(~S);
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
return EndomorphismBasis(GeoEndList, GensHf);

end intrinsic;
