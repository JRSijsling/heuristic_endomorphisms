/***
 *  Data per subfield
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function EndomorphismDataG3(GeoEndList, L, K);
// Input:   Algebraized, analytic, and rational representations of the basis of
//          an endomorphism ring of a Jacobian, along with its period matrix.
// Output:  Invariants for the associated ring, algebra, algebra tensored with
//          RR, and decompositions.

// The function always returns the endomorphism algebra; optionally, with
// AddTensor, AddRing and AddDecomposition, information is returned about the
// endomorphism ring (an order in the endomorphism algebra), the endomorphism
// algebra tensored with RR over QQ, and the idempotents needed for a
// decomposition of the isogeny class of the Jacobian.

// In the current version we use the compactified representation of these data,
// but the functions used can be applied in greater generality (this seemed like
// a good compromise).

// TODO: Make relative, and if possible work with general fields in the input so that other functions can go in the Wrapper.
        
if K eq L then
    // Geometric case:
    EndList0 := GeoEndList;
else
    Gp, Gf, Gphi := AutomorphismGroup(L);
    if Degree(K) eq 1 then
        GensHp := Generators(Gp);
    else
        Hp := FixedGroup(L, K);
        GensHp := Generators(Hp);
    end if;
    // Galois invariants and resulting endomorphisms:
    GensHf := [ Gphi(hp) : hp in GensHp ];
    EndList0 := EndomorphismBasisOverSubfield(GensHf, GeoEndList);
end if;
Rs := EndList0[3];

// Creation of relevant algebras:
g := #Rows(Rs[1]) div 2;
// Ambient matrix algebra, plus generators of the endomorphism ring:
A := Algebra(MatrixRing(Rationals(), 2*g)); GensA := [A ! Eltseq(R) : R in Rs];
// As a subalgebra:
B := sub<A | GensA>; GensB := [B ! gen : gen in GensA];
// As an associative algebra:
C := AssociativeAlgebra(B); GensC := [C ! gen : gen in GensB];

ED := [* *];
RepFactorsQQ := EndomorphismAlgebraFactorsQQG3(C);
Append(~ED, RepFactorsQQ);
RepFactorsRR := EndomorphismAlgebraFactorsRRG3(RepFactorsQQ);
Append(~ED, RepFactorsRR);
RepRing := EndomorphismRingG3(GensC, C);
Append(~ED, RepRing);
return ED;

end function;


function EndomorphismAlgebraFactorsQQG3(C);
// Determines the endomorphism algebra factors in the algorithm above. It
// heavily uses the classification. Also returns some other arguments needed for
// later use in determining the endomorphism ring and decompositions.

// Central decomposition:
Ds := DirectSumDecomposition(C);
RepFactorsQQ := [* *];
for D in Ds do
    RepFactorQQ := [* *];
    E1 := AlgebraOverCenter(D);
    F := IntegralRepresentationNF(BaseRing(E1));
    FDesc := Eltseq(MinimalPolynomial(F.1));
    E2 := ChangeRing(E1, F);
    test, d := IsSquare(Dimension(E2));
    if IsTotallyReal(F) then
        if d eq 1 then
            RepFactorQQ := [* "I", FDesc, d, -1 *];
        elif d eq 2 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Integers() ! Norm(DQFin);
            if not IsDefinite(Q) then
                RepFactorQQ := [* "II", FDesc, d, NDQ *];
            else
                // TODO: Case where we have a matrix ring over HH
                RepFactorQQ := [* "III", FDesc, d, NDQ *];
            end if;
        else
            RepFactorQQ := [* "II or III", FDesc, d, -1 *];
        end if;
    else
        if d eq 1 then
            RepFactorQQ := [* "IV", FDesc, d, -1 *];
        elif d eq 2 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Norm(DQFin);
            RepFactorQQ := [* "IV", FDesc, d, NDQ *];
        else
            RepFactorQQ := [* "IV", FDesc, d, -1 *];
        end if;
    end if;
    Append(~RepFactorsQQ, RepFactorQQ);
end for;

return RepFactorsQQ;

end function;


function EndomorphismAlgebraFactorsRRG3(RepFactorsQQ);
// Returns factors of endomorphism algebra tensored with RR.

RepFactorsRR := [];
for RepQQ in RepFactorsQQ do
    AlbertType := RepQQ[1];
    e := #RepQQ[2] - 1;
    if AlbertType eq "I" then
        RepFactorsRR cat:= [ "RR" : i in [1..e] ];
    elif AlbertType eq "II" then
        RepFactorsRR cat:= [ "M_2 (RR)" : i in [1..e] ];
    elif AlbertType eq "III" then
        RepFactorsRR cat:= [ "HH" : i in [1..e] ];
    elif AlbertType eq "OO or III" then
        RepFactorsRR cat:= [ "M2 (RR) or HH" : i in [1..e] ];
    elif AlbertType eq "IV" then
        d := RepQQ[3];
        if d eq 1 then
            RepFactorsRR cat:= [ "CC" : i in [1..(e div 2)] ];
        else
            str := Sprintf("M_%o (CC)", d);
            RepFactorsRR cat:= [ str : i in [1..(e div 2)] ];
        end if;
    end if;
end for;

return RepFactorsRR;

end function;


function EndomorphismRingG3(GensC, C);
// Returns endomorphism ring in the better ambient given by the previous
// algorithm.

OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM);
return Integers() ! ind;

end function;
