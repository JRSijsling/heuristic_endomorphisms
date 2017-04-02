/***
 *  Data per subfield
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


function EndomorphismData(GeoEndList, L, K, Gp, GensHp, Gphi, GeoFactorsQQ,
    Shorthand : AddTensor := false, AddRing := false, AddSatoTate := false,
    AddDecomposition := false);
// Input:   Algebraized, analytic, and rational representations of the basis of
//          an endomorphism ring of a Jacobian, along with its period matrix.
// Output:  Invariants for the associated ring, algebra, algebra tensored with
//          RR, and decompositions.

// The function always returns the endomorphism algebra; optionally, with
// AddTensor, AddRing and AddDecomposition, information is returned about the
// endomorphism ring (an order in the endomorphism algebra), the endomorphism
// algebra tensored with RR over QQ, and the idempotents needed for a decomposition
// of the isogeny class of the Jacobian.

// In the current version we use the compactified representation of these data,
// but the functions used can be applied in greater generality (this seemed like
// a good compromise).

// TODO: Make Gp, GensHp, Gphi keyword arguments and calculate them if not set.

if K eq L then
    // Geometric case:
    EndList0 := GeoEndList;
else
    // Galois invariants and resulting endomorphisms:
    GensHf := [ Gphi(hp) : hp in GensHp ];
    EndList0 := EndomorphismBasisOverSubfield(GensHf, GeoEndList);
end if;
Rs := EndList0[3];

// Creation of relevant algebras:
g := #Rows(Rs[1]) div 2;
// Ambient matrix algebra, plus generators of the endomorphism ring:
A := Algebra(MatrixRing(Rationals(), 2 * g));
GensA := [A ! Eltseq(R) : R in Rs];
// As a subalgebra:
B := sub<A | GensA>;
GensB := [B ! gen : gen in GensA];
// As an associative algebra:
C := AssociativeAlgebra(B);
GensC := [C ! gen : gen in GensB];

ED := [* *];
if g eq 2 then
    FactorsQQ, idemsC, CToPretty, RepFactorsQQ :=
        EndomorphismAlgebraFactorsQQG2(C);
    Append(~ED, RepFactorsQQ);
    if AddTensor then
        RepFactorsRR := EndomorphismAlgebraFactorsRRG2(FactorsQQ);
        Append(~ED, RepFactorsRR);
    end if;
    if AddRing then
        O, RepRing := EndomorphismRingG2(GensC, C, FactorsQQ, CToPretty);
        Append(~ED, RepRing);
    end if;
    if AddSatoTate then
        if K eq L then
            // Change parameters because we want this to work without Galois
            // group calculation.
            // TODO: Remove this hack: the fact that it works uses the
            // classification in an ugly way.
            RepSatoTate := SatoTateGroupG2(GeoEndRRShorthand(RepFactorsRR),
                RepFactorsRR, L, K, 0, 0, 0, GeoFactorsQQ, GeoEndList);
        else
            RepSatoTate := SatoTateGroupG2(Shorthand, RepFactorsRR, L, K, Gp,
                GensHp, Gphi, GeoFactorsQQ, GeoEndList);
        end if;
        Append(~ED, RepSatoTate);
    end if;
    if AddDecomposition then
        idems := DecompositionIdempotentsG2(idemsC, GensC, FactorsQQ, CToPretty,
            EndList0);
    else
        idems := [* [], [], [] *];
    end if;
    return ED, idems, FactorsQQ;
else
    // The generic non-genus 2 case:
    FactorsQQ := EndomorphismAlgebraFactorsGeneric(C);
    Append(~ED, FactorsQQ);
    if AddRing then
        // If asked for information about the order, we merely return its
        // discriminant:
        Append(~ED, EndomorphismRingGeneric(GensC));
    end if;
    // In general no quick option for AddTensor or AddDecomposition suggests itself,
    // so we do not give extra informations even if these flags are set.
    return ED, [], FactorsQQ;
end if;

end function;


function EndomorphismAlgebraFactorsGeneric(C);
// Determines endomorphism algebra factors in the generic case.
// Here we return the factors of the full endomorphism algebra, as algebras over
// their centers:

// Central decomposition:
Ds, idemsC := DirectSumDecomposition(C);
return [ AlgebraOverCenter(D) : D in Ds ];

end function;


function EndomorphismRingGeneric(GensC);
// Returns the endomorphism ring.
// TODO: This should be a better function that uses the reduced discriminant to
// return the index.

return Order(Integers(), GensC);

end function;


function EndomorphismAlgebraFactorsQQG2(C);
// Determines the endomorphism algebra factors in the algorithm above. It
// heavily uses the classification. Also returns some other arguments needed for
// later use in determining the endomorphism ring and decompositions.

// Central decomposition:
Ds, idemsC := DirectSumDecomposition(C);
issimple := #Ds eq 1;
if not issimple then
    // First case: central decomposition occurs.
    // Geometrically this means that we get non-isogenous elliptic curve
    // factors, so that the corresponding factors of the endomorphism algebra
    // are fields; we output our result as such.
    Es := [ AlgebraOverCenter(D) : D in Ds ];
    Fs := [* IntegralRepresentationNF(BaseRing(E)) : E in Es *];
    // FIXME: The -1 is for the uniformity of representation; I see no better
    // way.
    RepFactorsQQ := [ [* Eltseq(MinimalPolynomial(F.1)), -1 *] : F in Fs ];
    // The general procedure is to return the decomposition and describe its
    // simple elements. In this case there are two fields:
    return Fs, idemsC, 0, RepFactorsQQ;
else
    E1, f1 := AlgebraOverCenter(C);
    F, h := IntegralRepresentationNF(BaseRing(E1));
    // FIXME: As an aside, this functionality is weird since it seems to work
    // better without adding h. When changing ring it says that h is not a
    // homomorphism, which is freaky. Explanation?
    E2, f2 := ChangeRing(E1, F);
    if IsCommutative(E2) then
        // We first proceed to the subcase where the single factor is in fact a
        // field. Here we return information about the field as in the previous
        // case.
        RepFactorsQQ := [ [* Eltseq(MinimalPolynomial(F.1)), -1 *] ];
        return [ F ], idemsC, f1 * f2, RepFactorsQQ;
    else
        // Finally the case where the endomorphism ring is a quaternion algebra;
        // here we change its representation to be such:
        test, Q, f3 := IsQuaternionAlgebra(E2);
        DQFin := Discriminant(Q);
        NDQ := Norm(DQFin);
        RepFactorsQQ := [ [* Eltseq(MinimalPolynomial(F.1)), NDQ *] ];
        return [ Q ], idemsC, f1 * f2 * f3, RepFactorsQQ;
    end if;
end if;

end function;


function EndomorphismAlgebraFactorsRRG2(FactorsQQ);
// Returns factors of endomorphism algebra tensored with RR.

RepFactorsRR := [];
if #FactorsQQ eq 2 or IsCommutative(FactorsQQ[1]) then
    // Field case:
    for F in FactorsQQ do
        Inf := InfinitePlaces(F);
        for v in Inf do
            // FIXME: Guard against infinite places of QQ, which I do not know how
            // to treat uniformly:
            if Degree(F) eq 1 then
                Append(~RepFactorsRR, "RR");
            elif IsReal(v) then
                Append(~RepFactorsRR, "RR");
            else
                Append(~RepFactorsRR, "CC");
            end if;
        end for;
    end for;
    return RepFactorsRR;
else
    // Quaternion case:
    Q := FactorsQQ[1];
    // Used to be
    //DQFin, DQInf := Discriminant(Q);
    RamFin, RamInf := RamifiedPlaces(Q);
    F := BaseRing(Q);
    Inf := InfinitePlaces(F);
    for v in Inf do
        // FIXME: Guard against infinite places of QQ, which I do not know how
        // to treat uniformly:
        if Degree(F) eq 1 then
            if v in RamInf then
                Append(~RepFactorsRR, "HH");
            else
                Append(~RepFactorsRR, "M_2(RR)");
            end if;
        elif IsReal(v) then
            if v in RamInf then
                Append(~RepFactorsRR, "HH");
            else
                Append(~RepFactorsRR, "M_2(RR)");
            end if;
        else
            Append(~RepFactorsRR, "M_2(CC)");
        end if;
    end for;
    return RepFactorsRR;
end if;

end function;


function EndomorphismRingG2(GensC, C, FactorsQQ, CToPretty);
// Returns endomorphism ring in the better ambient given by the previous
// algorithm.

E1 := FactorsQQ[1];
F1 := BaseRing(E1);
if #FactorsQQ eq 2 or not IsCommutative(E1) then
    OC := Order(Integers(), GensC);
    DOC := Discriminant(OC);
    DOM := Discriminant(MaximalOrder(C));
    test, ind := IsSquare(DOC / DOM);
    RepRing := [ Integers() ! ind ];
    // In the case of a quaternion algebra over the rationals we can ask for a
    // bit more:
    if not IsCommutative(E1) and Degree(F1) eq 1 then
        // FIXME: Used to be
        //OE1 := Order([ CToPretty(gen) : gen in GensC ]);
        OE1 := QuaternionOrder([ CToPretty(gen) : gen in GensC ]);
        // FIXME: Next line deal with representation issues:
        if IsEichler(OE1) then
            Append(~RepRing, 1);
        else
            Append(~RepRing, 0);
        end if;
    else
        Append(~RepRing, -1);
    end if;
    return OC, RepRing;
else
    // Right now this functionality is only available for the case of a field,
    // since it does not work for quaternion algebras over general fields.
    // FIXME: Non-uniform approach for QQ:
    if Degree(E1) eq 1 then
        // Get ZZ
        return Integers(), [ 1, -1 ];
    end if;
    GensE1 := [ E1 ! CToPretty(gen) : gen in GensC ];
    OE1 := Order(GensE1);
    COE1 := Conductor(OE1);
    // Only return the norm for now:
    // TODO: Give factorization if it exists.
    RepRing := [ Norm(COE1), -1 ];
    //RepRing := [* Norm(COE1), [ Eltseq(E1 ! gen) : gen in Generators(COE1) ] *];
    return OE1, RepRing;
end if;

end function;


function DecompositionIdempotentsG2(idemsC, GensC, FactorsQQ, CToPretty,
    EndList0);
// Returns idempotents from which a decomposition can be reconstructed.
// The argument list is currently incredibly ugly because Magma refuses to
// return certain coercions as morphisms.

AsAlg0, As0, Rs0 := Explode(EndList0);
// Ignore the trivial central idempotent:
if #idemsC ne 2 then
    idems := [];
else
    idems := idemsC;
end if;
// First possibility: central decomposition of the algebra gives unique idempotents.
// If this occurs, then we already have these idempotents in GensC.
// Second possibility: obtained quaternion algebra isomorphic to matrix ring:
if not IsCommutative(FactorsQQ[1]) then
    Q := FactorsQQ[1];
    DQFin := Discriminant(Q);
    NDQ := Norm(DQFin);
    if NDQ eq 1 then
        PrettyToC := Inverse(CToPretty);
        GensQ := [ CToPretty(gen) : gen in GensC ];
        idems := [ PrettyToC(SmallIdempotent(Q, GensQ)) ];
    end if;
end if;
// Some ugly coercions are necessary to apply MatrixInBasis:
idems := [ [ Rationals() ! c : c in Eltseq(idem) ] : idem in idems ];
GensC := [ [ Rationals() ! c : c in Eltseq(gen) ] : gen in GensC ];
idemsRep := [ Eltseq(MatrixInBasis(idem, GensC)) : idem in idems ];
N := #AsAlg0;
idemsAlg := [ &+[ idem[n] * AsAlg0[n] : n in [1..N] ] : idem in idemsRep ];
idemsAn := [ &+[ idem[n] * As0[n] : n in [1..N] ] : idem in idemsRep ];
idemsRat := [ &+[ idem[n] * Rs0[n] : n in [1..N] ] : idem in idemsRep ];

return [* idemsAlg, idemsAn, idemsRat *];

end function;


function GeoEndRRShorthand(str);
// A bit of semantic sugar

case str:
    when ["RR"]: return "A";
    when ["RR", "RR"]: return "B";
    when ["RR", "CC"]: return "C";
    when ["CC", "RR"]: return "C";
    when ["CC", "CC"]: return "D";
    when ["M_2(RR)"]: return "E";
    when ["M_2(CC)"]: return "F";
    else: error Error("Shorthand algorithm obtains contradiction with classification");
end case;

end function;
