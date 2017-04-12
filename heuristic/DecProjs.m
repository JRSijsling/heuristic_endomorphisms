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


// TODO: Make possible without lattice?


intrinsic IdempotentDenominator(R::.) -> RngIntElt
{Degree of an idempotent or projection morphism.}

return LCM([ Denominator(c) : c in Eltseq(R) ]);

end intrinsic;


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
    // TODO: the usual nastyness with powers of quaternion algebra is again not
    // covered.
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

idems := [* idemsAlg, idemsRat, idemsAn *];
return idems;

end intrinsic;


intrinsic ProjectionFromIdempotent(idem::., P::.) -> List
{From an idempotent, extracts corresponding lattice and projection.}
    // Extract the complex field
    CC := BaseRing(P); RR := RealField(CC);

    // Create analytic idempotent and project:
    PEllHuge := P * idem;

    // Compute rank of projection
    gQuotient := NumericalRank(PEllHuge : Epsilon := RR`epsinv);

    PEllHugeSplit := SplitMatrix(PEllHuge);

    PEllBigSplit, col_numbers := SubmatrixOfRank(PEllHugeSplit, 2*gQuotient, "Rows"); // extract 2g independent rows
    PEllBigSplit := SaturateLattice(PEllHugeSplit, PEllBigSplit); // ensure that these rows generate the full lattice

    PEllBig := CombineMatrix(PEllBigSplit); // go back to the complex representation

    PreliminaryLatticeMatrix := SubmatrixOfRank(PEllBig, gQuotient, "Columns"); // extract g columns (i.e. decide which projection to use)

    PreliminaryLatticeMatrix := Transpose(PreliminaryLatticeMatrix); // necessary before calling SaturateLattice
    PEllBig := Transpose(PEllBig);
    LatticeMatrix := Transpose(SaturateLattice(PEllBig, PreliminaryLatticeMatrix));
end intrinsic;


intrinsic ProjectionFromIdempotents(idems::., P::.) -> List
{From idempotents, extracts corresponding lattices and projections.}

return [* ProjectionFromIdempotent(idem, P) : idem in idems *];

end intrinsic;
