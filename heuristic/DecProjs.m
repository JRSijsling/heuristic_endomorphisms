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


intrinsic IdempotentsFromLattice(Lat::List) -> .
{Finds idempotents over the smallest possible field by running over Lat.}

Lat := Reverse(Lat);
GeoEndoStruct := Lat[#Lat][2];
num_idemsgeo := NumberOfIdempotentsFromStructure(GeoEndoStruct);

n := #Lat; i := 1;
while i le n do
    K := Lat[i][1][2];
    EndoStruct := Lat[i][2];
    num_idems := NumberOfIdempotentsFromStructure(EndoStruct);
    if num_idems eq num_idemsgeo then
        break;
    end if;
    i +:= 1;
end while;

idems := IdempotentsFromStructure(EndoStruct);
return idems, K;

end intrinsic;


intrinsic NumberOfIdempotentsFromStructure(EndoStruct::List) -> RngIntElt
{Finds number of idempotents.}

EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
factors_QQ := EndoDesc[1];
num_idems := 0;
for factor_QQ in factors_QQ do
    albert, _, dim_sqrt, disc := Explode(factor_QQ);
    // TODO: the usual nastyness with powers of a quaternion algebra is again
    // not covered.
    if albert in ["II", "IV"] then
        if disc eq 1 then
            num_idems +:= dim_sqrt;
        end if;
    else
        num_idems +:= 1;
    end if;
end for;
return num_idems;

end intrinsic;


intrinsic IdempotentsFromStructure(EndoStruct::List) -> List
{Finds idempotents.}

g := #Rows(EndoStruct[1][1][1]);
EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
C, GensC := Explode(EndoAlg);
Ds := DirectSumDecomposition(C);
idemsC := &cat[ IdempotentsFromFactor(D, C, g) : D in Ds ];
idemsrep := MatricesFromIdempotents(idemsC, EndoStruct);
return idemsrep;

end intrinsic;


intrinsic IdempotentsFromFactor(D::., C::., RngIntElt) -> .
{Idempotents originating from the factor D.}

if g le 3 then
    return IdempotentsFromFactorG3(D, C);
else
    error "Higher genus not implemented yet";
end if;

end intrinsic;


intrinsic IdempotentsFromFactorG3(D::., C::.) -> .
{Idempotents originating from the factor D.}

// TODO: This 
E1, f1 := AlgebraOverCenter(D);
//F := ClearFieldDenominator(BaseRing(E1));
//if Type(F) eq FldNum then
//    F := OptimizedRepresentation(F);
//    F := ClearFieldDenominator(F);
//end if;
//E2, f2 := ChangeRing(E1, F);
E2 := E1;
if not IsCommutative(E2) then
    test_dim, d := IsSquare(Dimension(E2));
    if d eq 2 then
        test_quat, Q, f3 := IsQuaternionAlgebra(E2);
        if test_quat then
            test_mat, M, f4 := IsMatrixRing(Q : Isomorphism := true);
            //f := f1 * f2 * f3 * f4;
            f := f1 * f3 * f4;
            invf := Inverse(f);
            return [ C ! invf(M ! [1,0,0,0]), C ! invf(M ! [0,0,0,1]) ];
        end if;
    elif d eq 3 then
        idems_E2 := IdempotentsInMatrixAlgebra(E2);
        invf1 := Inverse(f1);
        return [ C ! invf1(idem_E2) : idem_E2 in idems_E2 ];
    else
        error "All cases in IdempotentsFromFactorG3 fell through";
    end if;
end if;
return [ C ! D ! 1 ];

end intrinsic;


intrinsic MatricesFromIdempotents(idems::SeqEnum, EndoStruct::List) -> SeqEnum
{Recovers matrices corresponding to idems.}

EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
GensTan := [ gen[1] : gen in EndoRep ];
GensHom := [ gen[2] : gen in EndoRep ];
GensApp := [ gen[3] : gen in EndoRep ];
C, GensC := Explode(EndoAlg);

idems := [ [ Rationals() ! c : c in Eltseq(idem) ] : idem in idems ];
GensC := [ [ Rationals() ! c : c in Eltseq(gen) ] : gen in GensC ];
idems := [ Eltseq(MatrixInBasis(idem, GensC)) : idem in idems ];

idemsRep := [ ];
for idem in idems do
    idemTan := &+[ idem[i] * GensTan[i] : i in [1..#GensTan] ];
    idemHom := &+[ idem[i] * GensHom[i] : i in [1..#GensHom] ];
    idemApp := &+[ idem[i] * GensApp[i] : i in [1..#GensApp] ];
    Append(~idemsRep, [* idemTan, idemHom, idemApp *]);
end for;
return idemsRep;

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
