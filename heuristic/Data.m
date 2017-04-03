/***
 *  Data per subfield
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

intrinsic EndomorphismData(EndList::SeqEnum) -> List
{Gives the endomorphism data over the subfield K, starting from.}

Rs := EndList[3];
// Creation of relevant algebras
g := #Rows(Rs[1]) div 2;
// Ambient matrix algebra, plus generators of the endomorphism ring
A := Algebra(MatrixRing(Rationals(), 2*g)); GensA := [ A ! Eltseq(R) : R in Rs ];
// As a subalgebra
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
// As an associative algebra
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

EndoDescList := [* *];
EndoDescQQ := EndomorphismDescriptionQQ(C);
Append(~EndoDescList, EndoDescQQ);
EndoDescRR := EndomorphismDescriptionRR(EndoDescQQ);
Append(~EndoDescList, EndoDescRR);
EndoDescZZ := EndomorphismDescriptionZZ(GensC, C);
Append(~EndoDescList, EndoDescZZ);

return C, EndoDescList;

end intrinsic;


intrinsic EndomorphismDescriptionQQ(C::AlgAss) -> .
{Decribes the factors of endomorphism algebra.}

// Central decomposition
Ds := DirectSumDecomposition(C);
EndoDescQQ := [* *];
for D in Ds do
    DescFactorQQ := [* *];
    E1 := AlgebraOverCenter(D);
    F := IntegralRepresentationNF(BaseRing(E1));
    FDesc := Eltseq(MinimalPolynomial(F.1));
    E2 := ChangeRing(E1, F);
    test, d := IsSquare(Dimension(E2));
    if IsTotallyReal(F) then
        if d eq 1 then
            DescFactorQQ := [* "I", FDesc, d, -1 *];
        elif d eq 2 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Integers() ! Norm(DQFin);
            if not IsDefinite(Q) then
                DescFactorQQ := [* "II", FDesc, d, NDQ *];
            else
                // TODO: Case where we have a matrix ring over HH
                DescFactorQQ := [* "III", FDesc, d, NDQ *];
            end if;
        else
            DescFactorQQ := [* "II or III", FDesc, d, -1 *];
        end if;
    else
        if d eq 1 then
            DescFactorQQ := [* "IV", FDesc, d, -1 *];
        elif d eq 2 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Norm(DQFin);
            DescFactorQQ := [* "IV", FDesc, d, NDQ *];
        else
            DescFactorQQ := [* "IV", FDesc, d, -1 *];
        end if;
    end if;
    Append(~EndoDescQQ, DescFactorQQ);
end for;
return EndoDescQQ;

end intrinsic;


intrinsic EndomorphismDescriptionRR(EndoDescQQ::List) -> .
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
        EndoDescRR cat:= [ "M2 (RR) or HH" : i in [1..e] ];
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
return EndoDescRR;

end intrinsic;


intrinsic EndomorphismDescriptionZZ(GensC::SeqEnum, C::AlgAss) -> .
{Describes of the endomorphism ring.}

OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM);
return Integers() ! ind;

end intrinsic;
