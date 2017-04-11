/***
 *  Subrings of the geometric endomorphism ring over subfields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic EndomorphismBasis(GeoEndList::List, GalK::List) -> List
{Extracts basis over subfield determines by a list of automorphisms.}

AsAlg, Rs, As := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);
GensH, Gphi := Explode(GalK);
K := FixedField(L, [ Gphi(gen) : gen in GensH ]);

if #GensH eq 0 then
    return GeoEndList;
end if;

// The vector space representing the full endomorphism algebra
n := #AsAlg;
Ker := VectorSpace(Rationals(), n);
// Successively filter by writing down the conditions for a matrix to be fixed
// under a given generator
for gen in GensH do
    sigma := Gphi(gen);
    Msigma := Matrix([MatrixInBasis(ConjugateMatrix(sigma, A), AsAlg) : A in AsAlg]);
    Msigma -:= IdentityMatrix(Rationals(), n);
    Ker meet:= Kernel(Msigma);
end for;

// Retrieve representations for a basis of of the endomorphism ring by taking a
// pure lattice / saturation
Lat := PureLattice(Lattice(Matrix(Basis(Ker))));
B := Basis(Lat);

// Constructing said basis
AsAlg := [ &+[ b[i] * AsAlg[i] : i in [1..n] ] : b in B ];
Rs    := [ &+[ b[i] * Rs[i]    : i in [1..n] ] : b in B ];
As    := [ &+[ b[i] * As[i]    : i in [1..n] ] : b in B ];
AsAlg := [ Matrix(K, A) : A in AsAlg ];
return [* AsAlg, Rs, As *];

end intrinsic;


intrinsic EndomorphismBasis(GeoEndList::List, K::Fld) -> List
{Extracts basis over a general field.}

AsAlg, As, Rs := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);

GalK := SubgroupGeneratorsUpToConjugacy(L, K);
return EndomorphismBasis(GeoEndList, GalK);

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


intrinsic SubgroupGeneratorsUpToConjugacy(L::Fld, K::Fld) -> List
{Finds the subgroup generators up to conjugacy that correspond to the intersection of L and K, where L is Galois.}

if (Degree(K) eq 1) or (Type(L) eq FldRat) then
    Gp, Gf, Gphi := AutomorphismGroup(L);
    return [* Generators(Gp), Gphi *];
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
    GalK := [* Generators(FixedGroup(L, KL)), Gphi *];
else
    // Case where we need to go down all the way
    GalK := [* Generators(Gp), Gphi *];
end if;

return GalK;

end intrinsic;
