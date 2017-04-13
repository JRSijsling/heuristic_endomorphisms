/***
 *  Elliptic curves from decompositions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


// TODO: Make possible without lattice? Transpose? Latter prolly not.


intrinsic IdempotentsFromLattice(Lat::List) -> .
{Finds idempotents over the smallest possible field by running over Lat.}

LatReps, LatAlgs, LatDescs := Explode(Lat);
GeoLatRep := LatReps[#LatReps]; GeoLatAlg := LatAlgs[#LatReps]; GeoLatDesc := LatDescs[#LatReps];
num_idems_geo := NumberOfIdempotentsFromLatticeEntry(GeoLatRep, GeoLatAlg, GeoLatDesc);

// Potentially computes final entry two times but who cares
n := #LatReps; i := 1;
while i le n do
    LatRep := LatReps[i]; LatAlg := LatAlgs[i]; LatDesc := LatDescs[i];
    num_idems := NumberOfIdempotentsFromLatticeEntry(LatRep, LatAlg, LatDesc);
    if num_idems eq num_idems_geo then
        break;
    end if;
    i +:= 1;
end while;

idems := IdempotentsFromLatticeEntry(LatRep, LatAlg, LatDesc);
K := BaseRing(idems[1][1]);
return idems, K;

end intrinsic;


intrinsic NumberOfIdempotentsFromLatticeEntry(LatRep::List, LatAlg::List, LatDesc::List) -> List
{Finds number of idempotents.}

K, factors_QQ := Explode(LatDesc);
num_idems := 0;
for factor_QQ in factors_QQ do
    albert, _, d, disc := Explode(factor_QQ);
    // TODO: the usual nastyness with powers of a quaternion algebra is again not covered.
    if albert in ["II", "IV"] then
        if disc eq 1 then
            num_idems +:= d;
        end if;
    else
        num_idems +:= 1;
    end if;
end for;
return num_idems;

end intrinsic;


intrinsic IdempotentsFromLatticeEntry(LatRep::List, LatAlg::List, LatDesc::List) -> List
{Finds idempotents.}

K, AsAlg, Rs, As := Explode(LatRep);
K, C, GensC := Explode(LatAlg);

Ds := DirectSumDecomposition(C);
idems_C := &cat[ IdempotentsFromFactor(D, C) : D in Ds ];
idems_rep := MatricesFromIdempotents(idems_C, LatRep, LatAlg, LatDesc);
return idems_rep;

end intrinsic;


intrinsic IdempotentsFromFactor(D::., C::.) -> .
{Idempotents originating from the factor D.}

// TODO: Larger d (this is actually an issue in genus 3)
E1, f1 := AlgebraOverCenter(D);
F := ClearFieldDenominator(BaseRing(E1));
if Type(F) eq FldNum then
    F := OptimizedRepresentation(F);
    F := ClearFieldDenominator(F);
end if;
E2, f2 := ChangeRing(E1, F);
test_dim, d := IsSquare(Dimension(E2));
if test_dim then
    if d eq 2 then
        test_quat, Q, f3 := IsQuaternionAlgebra(E2);
        if test_quat then
            test_mat, M, f4 := IsMatrixRing(Q : Isomorphism := true);
            f := f1 * f2 * f3 * f4;
            invf := Inverse(f);
            return [ C ! invf(M![1,0,0,0]), C ! invf(M![0,0,0,1]) ];
        end if;
    end if;
end if;
return [ C ! D ! 1 ];

end intrinsic;


intrinsic MatricesFromIdempotents(idems, LatRep, LatAlg, LatDesc) -> List
{Recovers matrices corresponding to idems.}

K, AsAlg, Rs, As := Explode(LatRep);
K, C, GensC := Explode(LatAlg);

idemsC := [ [ Rationals() ! c : c in Eltseq(idem) ] : idem in idems ];
GensC := [ [ Rationals() ! c : c in Eltseq(gen) ] : gen in GensC ];
idemsRep := [ Eltseq(MatrixInBasis(idemC, GensC)) : idemC in idemsC ];

n := #AsAlg;
idemsAlg := [ &+[ idemRep[i] * AsAlg[i] : i in [1..n] ] : idemRep in idemsRep ];
idemsRat := [ &+[ idemRep[i] * Rs[i] : i in [1..n] ] : idemRep in idemsRep ];
idemsAn := [ &+[ idemRep[i] * As[i] : i in [1..n] ] : idemRep in idemsRep ];
return [* idemsAlg, idemsRat, idemsAn *];

end intrinsic;


intrinsic ProjectionFromIdempotent(P::., idem::List) -> List
{From an idempotent, extracts corresponding lattice and projection.}
    
idem_alg, idem_rat, idem_app := Explode(idem);
// Extract the complex field
CC := BaseRing(P); RR := RealField(CC);

// Create analytic idempotent and project:
PEllHuge := P * idem_app;

// Compute rank of projection
gQuotient := NumericalRank(PEllHuge : Epsilon := RR`epsinv);

PEllHugeSplit := SplitMatrix(PEllHuge);

PEllBigSplit, col_numbers := SubmatrixOfRank(PEllHugeSplit, 2*gQuotient : ColumnsOrRows := "Rows"); // extract 2g independent rows
PEllBigSplit := SaturateLattice(PEllHugeSplit, PEllBigSplit); // ensure that these rows generate the full lattice

PEllBig := CombineMatrix(PEllBigSplit, CC); // go back to the complex representation

PreliminaryLatticeMatrix, s0 := SubmatrixOfRank(PEllBig, gQuotient : ColumnsOrRows := "Columns"); // extract g columns (i.e. decide which projection to use)

PreliminaryLatticeMatrix := Transpose(PreliminaryLatticeMatrix); // necessary before calling SaturateLattice
PEllBig := Transpose(PEllBig);
LatticeMatrix := Transpose(SaturateLattice(PEllBig, PreliminaryLatticeMatrix));

rows_alg := Rows(idem_alg); proj_alg := Matrix([ rows_alg[i] : i in s0 ]);
rows_rat := Rows(idem_rat); proj_rat := Matrix([ rows_rat[i] : i in s0 ]);
rows_app := Rows(idem_app); proj_app := Matrix([ rows_app[i] : i in s0 ]);
proj := [* proj_alg, proj_rat, proj_app *];

return [* LatticeMatrix, proj *];
end intrinsic;


intrinsic ProjectionsFromIdempotents(P::., idems::.) -> List
{From idempotents, extracts corresponding lattices and projections.}

lats_projs := [* *];
for i:=1 to #idems[1] do
    idem_alg := idems[1][i]; idem_rat := idems[2][i]; idem_app := idems[3][i];
    idem := [* idem_alg, idem_rat, idem_app *];
    Append(~lats_projs, ProjectionFromIdempotent(P, idem));
end for;
return lats_projs;

end intrinsic;


intrinsic IdempotentDenominator(R::.) -> RngIntElt
{Degree of an idempotent or projection morphism.}

return LCM([ Denominator(c) : c in Eltseq(R) ]);

end intrinsic;
