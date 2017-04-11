/***
 *  Structural description per subfield
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic EndomorphismStructure(GeoEndList::List, GensHf::SeqEnum) -> List
{Gives the endomorphism structure over the subfield corresponding to GensHf, starting from a list of representations.}

EndoReps := EndomorphismBasis(GeoEndList, GensHf);
EndoStruct, EndoDesc := EndomorphismStructureAndDescription(EndoReps);
SatoTate := SatoTateGroup(GeoEndList, GensHf);
Append(~EndoStruct, SatoTate); Append(~EndoDesc, SatoTate);
return EndoReps, EndoStruct, EndoDesc;

end intrinsic;


intrinsic EndomorphismStructure(GeoEndList::List, K::Fld) -> List
{Gives the endomorphism structure over the subfield K, starting from a list of representations.}

EndoReps := EndomorphismBasis(GeoEndList, K);
EndoStruct, EndoDesc := EndomorphismStructureAndDescription(EndoReps);
SatoTate := SatoTateGroup(GeoEndList, K);
Append(~EndoStruct, SatoTate); Append(~EndoDesc, SatoTate);
return EndoReps, EndoStruct, EndoDesc;

end intrinsic;


intrinsic EndomorphismStructureAndDescription(EndList::List) -> List
{Gives the endomorphism structure, starting from a list of representations.}

Rs := EndList[2];
// Creation of relevant algebras
g := #Rows(Rs[1]) div 2;
// Ambient matrix algebra, plus generators of the endomorphism ring
A := Algebra(MatrixRing(Rationals(), 2*g)); GensA := [ A ! Eltseq(R) : R in Rs ];
// As a subalgebra
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
// As an associative algebra
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

EndoStruct := [* *]; EndoDesc := [* *];
EndoStructQQ, EndoDescQQ := EndomorphismStructureQQ(C);
Append(~EndoStruct, EndoStructQQ); Append(~EndoDesc, EndoDescQQ);
EndoStructZZ, EndoDescZZ := EndomorphismStructureZZ(C, GensC);
Append(~EndoStruct, EndoStructZZ); Append(~EndoDesc, EndoDescZZ);
EndoStructRR, EndoDescRR := EndomorphismStructureRR(C, EndoDescQQ);
Append(~EndoStruct, EndoStructRR); Append(~EndoDesc, EndoDescRR);

return EndoStruct, EndoDesc;

end intrinsic;


intrinsic EndomorphismStructureQQ(C::AlgAss : Optimize := true) -> .
{Decribes the factors of endomorphism algebra.}
// FIXME: Right now a non-Optimized version is not yet easily accessible from
// Sage. On the other hand, it is not a big operation, so I hav not spent time
// on this.

// Central decomposition
Ds := DirectSumDecomposition(C);
EndoDescQQ := [* *];
for D in Ds do
    DescFactorQQ := [* *];
    E1 := AlgebraOverCenter(D);
    F := ClearFieldDenominator(BaseRing(E1));
    if (Type(F) eq FldNum and Optimize) then
        F := OptimizedRepresentation(F);
    end if;
    FDesc := Eltseq(MinimalPolynomial(F.1));
    // Next line is not strictly necessary, but for uniformity of the description:
    FDesc := [ Integers() ! c : c in FDesc ];
    E2 := ChangeRing(E1, F);
    test, d := IsSquare(Dimension(E2));
    if IsTotallyReal(F) then
        if d eq 1 then
            DescFactorQQ := [* "I", FDesc, d, 1 *];
        elif d eq 2 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Integers() ! Norm(DQFin);
            if not IsDefinite(Q) then
                DescFactorQQ := [* "II", FDesc, d, NDQ *];
            else
                DescFactorQQ := [* "III", FDesc, d, NDQ *];
            end if;
        elif d eq 3 then
            DescFactorQQ := [* "II", FDesc, d, -1 *];
        else
            /* FIXME: We do not know what happens here, even when using the
             * extended Albert classification that I have applied. Testing for
             * a matrix ring can be done with
             *     Norm(Discriminant(MaximalOrder(E2)));
             * but that is only a necessary condition; there may be
             * ramification at infinity only, in which case this does not tell
             * enough. For now we get by; this should be addressed with more
             * general functionality for algebras, not by our package. */
            DescFactorQQ := [* "II or III", FDesc, d, -1 *];
        end if;
    else
        if d eq 1 then
            DescFactorQQ := [* "IV", FDesc, d, 1 *];
        elif d eq 2 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Norm(DQFin);
            DescFactorQQ := [* "IV", FDesc, d, NDQ *];
        elif d eq 3 then
            DescFactorQQ := [* "IV", FDesc, d, 1 *];
        else
            DescFactorQQ := [* "IV", FDesc, d, -1 *];
        end if;
    end if;
    Append(~EndoDescQQ, DescFactorQQ);
end for;
return C, EndoDescQQ;

end intrinsic;


intrinsic EndomorphismStructureRR(C::AlgAss, EndoDescQQ::List) -> .
{Decribes the factors of endomorphism algebra tensored with RR.}

EndoDescRR := [ ];
for DescFactorQQ in EndoDescQQ do
    AlbertType := DescFactorQQ[1];
    e := #DescFactorQQ[2] - 1;
    if AlbertType eq "I" then
        EndoDescRR cat:= [ "RR" : i in [1..e] ];
    elif AlbertType eq "II" then
        EndoDescRR cat:= [ "M_2 (RR)" : i in [1..e] ];
    elif AlbertType eq "III" then
        EndoDescRR cat:= [ "HH" : i in [1..e] ];
    elif AlbertType eq "OO or III" then
        EndoDescRR cat:= [ "M_2 (RR) or HH" : i in [1..e] ];
    elif AlbertType eq "IV" then
        d := DescFactorQQ[3];
        if d eq 1 then
            EndoDescRR cat:= [ "CC" : i in [1..(e div 2)] ];
        else
            str := Sprintf("M_%o (CC)", d);
            EndoDescRR cat:= [ str : i in [1..(e div 2)] ];
        end if;
    end if;
end for;
return EndoDescRR, EndoDescRR;

end intrinsic;


intrinsic EndomorphismStructureZZ(C::AlgAss, GensC::SeqEnum : Optimize := true) -> .
{Describes of the endomorphism ring.}

// Calculating index
OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM);

// Test whether Eichler in a quaternion algebra
Ds := DirectSumDecomposition(C);
if #Ds eq 1 then
    E1, f1 := AlgebraOverCenter(C);
    F := ClearFieldDenominator(BaseRing(E1));
    if (Type(F) eq FldNum and Optimize) then
        F := OptimizedRepresentation(F);
    end if;
    E2, f2 := ChangeRing(E1, F);
    test, d := IsSquare(Dimension(E2));
    if d eq 2 then
        test, Q, f3 := IsQuaternionAlgebra(E2);
        if test then
            f := f1 * f2 * f3;
            OO := QuaternionOrder([ f(gen) : gen in GensC ]);
            if IsEichler(OO) then
                return OC, [ Integers() ! ind, 1 ];
            end if;
        end if;
    end if;
end if;
return OC, [ Integers() ! ind, -1 ];

end intrinsic;
