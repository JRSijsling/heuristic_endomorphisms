/***
 *  Determining the Sato-Tate group
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic SatoTateGroupG2(GeoEndList::List, GalK::List) -> MonStgElt
{Determines the Sato-Tate group in genus 2.}

// The first two steps are potentially inefficient, but in practice not that
// bad, so I keep it in.
// Alternatively, we can do this by keeping track of keyword arguments:
// (1) Pass along desc_RR from a previous calculation of EndoAlg, which
//     typically comes just one step before;
// (2) Pass along desc_RR_geo when using lattices by first calculating the
//     geometric case and using the result.
// The case where everything is defined over the base field will then still
// calculate things twice. This is hopefully forgiveable; alternatively, we can
// cover it by some trivial return statements in the EndomorphismStructure
// functionality.

GeoEndoAlg, GeoEndoDesc := EndomorphismAlgebraAndDescription(GeoEndList);
desc_RR_geo := GeoEndoAlg[3];
Shorthand := SatoTateShorthandG2(desc_RR_geo);

EndoReps := EndomorphismBasis(GeoEndList, GalK);
EndoAlg, EndoDesc := EndomorphismAlgebraAndDescription(EndoReps);
desc_RR := EndoAlg[3];

AsAlg, Rs, As := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);
GensH, Gphi := Explode(GalK);
H := sub< Domain(Gphi) | GensH >;
K := FixedField(L, H);

if Shorthand eq "A" then
    return "USp(4)";

elif Shorthand eq "B" then
    if desc_RR eq ["RR"] then
        return "N(G_{3,3})";
    elif desc_RR eq ["RR", "RR"] then
        return "G_{3,3}";
    end if;

elif Shorthand eq "C" then
    if desc_RR eq ["RR", "RR"] then
        return "N(G_{1,3})";
    elif desc_RR eq ["RR", "CC"] or desc_RR eq ["CC", "RR"] then
        return "G_{1,3}";
    end if;

elif Shorthand eq "D" then
    if desc_RR eq ["RR"] then
        return "F_{ac}";
    elif desc_RR eq ["RR", "RR"] then
        if IsIsomorphic(H, CyclicGroup(2)) then
            return "F_{ab}";
        elif IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "F_{a,b}";
        end if;
    elif desc_RR eq ["RR", "CC"] or desc_RR eq ["CC", "RR"] then
        return "F_a";
    elif desc_RR eq ["CC", "CC"] then
        return "F";
    end if;

elif Shorthand eq "E" then
    if desc_RR eq ["RR"] then
        // Magma being silly with DihedralGroup(2)...
        if IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "J(E_2)";
        elif IsIsomorphic(H, DihedralGroup(3)) then
            return "J(E_3)";
        elif IsIsomorphic(H, DihedralGroup(4)) then
            return "J(E_4)";
        elif IsIsomorphic(H, DihedralGroup(6)) then
            return "J(E_6)";
        end if;
    elif desc_RR eq ["RR", "RR"] then
        return "J(E_1)";
    elif desc_RR eq ["CC"] then
        if IsIsomorphic(H, CyclicGroup(2)) then
            return "E_2";
        elif IsIsomorphic(H, CyclicGroup(3)) then
            return "E_3";
        elif IsIsomorphic(H, CyclicGroup(4)) then
            return "E_4";
        elif IsIsomorphic(H, CyclicGroup(6)) then
            return "E_6";
        end if;
    elif desc_RR eq ["M_2 (RR)"] then
        return "E_1";
    end if;

elif Shorthand eq "F" then
    if desc_RR eq ["RR"] then
        if IsIsomorphic(H, DihedralGroup(4)) then
            return "D_{4,1}";

        elif IsIsomorphic(H, DihedralGroup(6)) then
            // See Fité--Kedlaya--Rotger--Sutherland (4.3) for the next step
            H_prime := Center(H); GensH_prime := Generators(H_prime);
            GalK_prime := [* GensH_prime, Gphi *];
            EndoReps_prime := EndomorphismBasis(GeoEndList, GalK_prime);
            EndoAlg_prime, EndoDesc_prime := EndomorphismAlgebraAndDescription(EndoReps_prime);
            desc_RR_prime := EndoAlg_prime[3];
            if desc_RR_prime eq ["M_2 (RR)"] then
                return "D_{6,1}";
            elif desc_RR_prime eq ["HH"] then
                return "J(D_3)";
            end if;

        elif IsIsomorphic(H, SymmetricGroup(4)) then
            return "O_1";
        elif IsIsomorphic(H, DirectProduct(DirectProduct(CyclicGroup(2), CyclicGroup(2)), CyclicGroup(2))) then
            return "J(D_2)";
        elif IsIsomorphic(H, DirectProduct(DihedralGroup(4), CyclicGroup(2))) then
            return "J(D_4)";
        elif IsIsomorphic(H, DirectProduct(DihedralGroup(6), CyclicGroup(2))) then
            return "J(D_6)";
        elif IsIsomorphic(H, DirectProduct(AlternatingGroup(4), CyclicGroup(2))) then
            return "J(T)";
        elif IsIsomorphic(H, DirectProduct(SymmetricGroup(4), CyclicGroup(2))) then
            return "J(O)";
        end if;

    elif desc_RR eq ["RR", "RR"] then
        if IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "D_{2,1}";
        elif IsIsomorphic(H, DihedralGroup(3)) then
            return "D_{3,2}";
        elif IsIsomorphic(H, DihedralGroup(4)) then
            return "D_{4,2}";
        elif IsIsomorphic(H, DihedralGroup(6)) then
            return "D_{6,2}";
        end if;

    elif desc_RR eq ["CC"] then
        if IsIsomorphic(H, CyclicGroup(4)) then
            return "C_{4,1}";

        elif IsIsomorphic(H, CyclicGroup(6)) then
            // Here we take the unique subgroup of order 2:
            H_prime := Subgroups(H : OrderEqual := 2)[1]`subgroup; GensH_prime := Generators(H_prime);
            GalK_prime := [* GensH_prime, Gphi *];
            EndoReps_prime := EndomorphismBasis(GeoEndList, GalK_prime);
            EndoAlg_prime, EndoDesc_prime := EndomorphismAlgebraAndDescription(EndoReps_prime);
            desc_RR_prime := EndoAlg_prime[3];
            if desc_RR_prime eq ["M_2 (RR)"] then
                return "C_{6,1}";
            elif desc_RR_prime eq ["HH"] then
                return "J(C_3)";
            end if;

        elif IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            // In this case it suffices to check whether the polynomial that
            // defines the center of the geometric endomorphism ring in fact
            // has a root in the ground field:
            f := DefiningPolynomial(BaseRing(GeoEndoAlg[1]));
            if HasRoot(f, K) then
                return "D_2";
            else
                return "J(C_2)";
            end if;

        elif IsIsomorphic(H, DirectProduct(CyclicGroup(4), CyclicGroup(2))) then
            return "J(C_4)";
        elif IsIsomorphic(H, DirectProduct(CyclicGroup(6), CyclicGroup(2))) then
            return "J(C_6)";
        elif IsIsomorphic(H, DihedralGroup(3)) then
            return "D_3";
        elif IsIsomorphic(H, DihedralGroup(4)) then
            return "D_4";
        elif IsIsomorphic(H, DihedralGroup(6)) then
            return "D_6";
        elif IsIsomorphic(H, AlternatingGroup(4)) then
            return "T";
        elif IsIsomorphic(H, SymmetricGroup(4)) then
            return "O";
        end if;

    elif desc_RR eq ["CC", "CC"] then
        if IsIsomorphic(H, CyclicGroup(2)) then
            return "C_2";
        elif IsIsomorphic(H, CyclicGroup(3)) then
            return "C_3";
        elif IsIsomorphic(H, CyclicGroup(4)) then
            return "C_4";
        elif IsIsomorphic(H, CyclicGroup(6)) then
            return "C_6";
        end if;

    elif desc_RR eq ["M_2 (RR)"] then
        return "C_{2,1}";
    elif desc_RR eq ["HH"] then
        return "J(C_1)";
    elif desc_RR eq ["M_2 (CC)"] then
        return "C_1";
    end if;
end if;
error Error("All cases in SatoTateGroupG2 fell through");

end intrinsic;


intrinsic SatoTateShorthandG2(desc_RR::SeqEnum) -> MonStgElt
{Finds the letter describing the neutral connected component of the Sato-Tate group.}

case desc_RR:
    when ["RR"]:       return "A";
    when ["RR", "RR"]: return "B";
    when ["RR", "CC"]: return "C";
    when ["CC", "RR"]: return "C";
    when ["CC", "CC"]: return "D";
    when ["M_2 (RR)"]: return "E";
    when ["M_2 (CC)"]: return "F";
    else: error Error("Shorthand algorithm obtains contradiction with classification");
end case;

end intrinsic;


intrinsic SatoTateGroup(GeoEndList::List, GalK::List) -> List
{Determines Sato-Tate group by subgroup of Galois group.}

AsAlg, Rs, As := Explode(GeoEndList);
g := #Rows(AsAlg[1]);

if g eq 2 then
    return SatoTateGroupG2(GeoEndList, GalK);
else
    // TODO: Add other cases when they appear.
    return "";
end if;

end intrinsic;


intrinsic SatoTateGroup(GeoEndList::List, K::Fld) -> List
{Determines Sato-Tate group over K.}

AsAlg, As, Rs := Explode(GeoEndList);
L := BaseRing(AsAlg[1]);

GalK := SubgroupGeneratorsUpToConjugacy(L, K);
return SatoTateGroup(GeoEndList, GalK);

end intrinsic;
