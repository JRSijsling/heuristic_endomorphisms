/***
 *  Determining the Sato-Tate group
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

function SatoTateGroupG2(Shorthand, FactorsRR, L, K, Gp, GensHp, Gphi,
    GeoFactorsQQ, GeoEndList);
// Returns the Sato-Tate group by an enormous case distinction. The amount of
// variables that this needs is outrageous, but I cannot do better for now.
// TODO: Make Gp, Hp, Gphi keyword arguments, calculate if not set.

if Shorthand eq "A" then
    return "USp(4)";
elif Shorthand eq "B" then
    if FactorsRR eq ["RR"] then
        return "N(G_{3,3})";
    elif FactorsRR eq ["RR", "RR"] then
        return "G_{3,3}";
    end if;
elif Shorthand eq "C" then
    if FactorsRR eq ["RR", "RR"] then
        return "N(G_{1,3})";
    elif FactorsRR eq ["RR", "CC"] or FactorsRR eq ["CC", "RR"] then
        return "G_{1,3}";
    end if;
elif Shorthand eq "D" then
    if FactorsRR eq ["RR"] then
        return "F_{ac}";
    elif FactorsRR eq ["RR", "RR"] then
        Hp := sub<Gp | GensHp>;
        if IsIsomorphic(Hp, CyclicGroup(2)) then
            return "F_{ab}";
        elif IsIsomorphic(Hp, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "F_{a,b}";
        end if;
    elif FactorsRR eq ["RR", "CC"] or FactorsRR eq ["CC", "RR"] then
        return "F_a";
    elif FactorsRR eq ["CC", "CC"] then
        return "F";
    end if;
elif Shorthand eq "E" then
    if FactorsRR eq ["RR"] then
        Hp := sub<Gp | GensHp>;
        // Magma being silly with DihedralGroup(2)...
        if IsIsomorphic(Hp, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "J(E_2)";
        elif IsIsomorphic(Hp, DihedralGroup(3)) then
            return "J(E_3)";
        elif IsIsomorphic(Hp, DihedralGroup(4)) then
            return "J(E_4)";
        elif IsIsomorphic(Hp, DihedralGroup(6)) then
            return "J(E_6)";
        end if;
    elif FactorsRR eq ["RR", "RR"] then
        return "J(E_1)";
    elif FactorsRR eq ["CC"] then
        Hp := sub<Gp | GensHp>;
        if IsIsomorphic(Hp, CyclicGroup(2)) then
            return "E_2";
        elif IsIsomorphic(Hp, CyclicGroup(3)) then
            return "E_3";
        elif IsIsomorphic(Hp, CyclicGroup(4)) then
            return "E_4";
        elif IsIsomorphic(Hp, CyclicGroup(6)) then
            return "E_6";
        end if;
    elif FactorsRR eq ["M_2(RR)"] then
        return "E_1";
    end if;
elif Shorthand eq "F" then
    if FactorsRR eq ["RR"] then
        Hp := sub<Gp | GensHp>;
        if IsIsomorphic(Hp, DihedralGroup(4)) then
            return "D_{4,1}";
        elif IsIsomorphic(Hp, DihedralGroup(6)) then
            // We split the group as in Fité--Kedlaya--Rotger--Sutherland (4.3).
            // It suffices to look at the subgroups of order 2 and find one for
            // which the polynomial defining the center of the geometric
            // endomorphism ring does not have a root in the corresponding
            // field:
            SubHp := Center(Hp);
            SubHf := [ Gphi(g) : g in SubHp ];
            Kint := FixedField(L, SubHf);
            GensHintp := Generators(SubHp);
            ED := EndomorphismData(GeoEndList, L, Kint, Gp, GensHintp, Gphi,
                GeoFactorsQQ, Shorthand : AddTensor := true, AddRing := false,
                AddSatoTate := false, AddDecomposition := false);
            if ED[2] eq ["M_2(RR)"] then
                return "D_{6,1}";
            elif ED[2] eq ["HH"] then
                return "J(D_3)";
            end if;
        elif IsIsomorphic(Hp, SymmetricGroup(4)) then
            return "O_1";
        elif IsIsomorphic(Hp, DirectProduct(DirectProduct(CyclicGroup(2), CyclicGroup(2)), CyclicGroup(2))) then
            return "J(D_2)";
        elif IsIsomorphic(Hp, DirectProduct(DihedralGroup(4), CyclicGroup(2))) then
            return "J(D_4)";
        elif IsIsomorphic(Hp, DirectProduct(DihedralGroup(6), CyclicGroup(2))) then
            return "J(D_6)";
        elif IsIsomorphic(Hp, DirectProduct(AlternatingGroup(4), CyclicGroup(2))) then
            return "J(T)";
        elif IsIsomorphic(Hp, DirectProduct(SymmetricGroup(4), CyclicGroup(2))) then
            return "J(O)";
        end if;
    elif FactorsRR eq ["RR", "RR"] then
        Hp := sub<Gp | GensHp>;
        if IsIsomorphic(Hp, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "D_{2,1}";
        elif IsIsomorphic(Hp, DihedralGroup(3)) then
            return "D_{3,2}";
        elif IsIsomorphic(Hp, DihedralGroup(4)) then
            return "D_{4,2}";
        elif IsIsomorphic(Hp, DihedralGroup(6)) then
            return "D_{6,2}";
        end if;
    elif FactorsRR eq ["CC"] then
        Hp := sub<Gp | GensHp>;
        if IsIsomorphic(Hp, CyclicGroup(4)) then
            return "C_{4,1}";
        elif IsIsomorphic(Hp, CyclicGroup(6)) then
            // This case is simpler, as there is only one subgroup of order 2:
            Hintp := Subgroups(Hp : OrderEqual := 2)[1]`subgroup;
            GensHintp := Generators(Hintp);
            GensHintf := [ Gphi(g) : g in GensHintp ];
            Kint := FixedField(L, GensHintf);
            if Degree(L) ne Degree(Kint) * #Hintp then
                // FIXME: Further bug guard:
                error Error("Magma took an incorrect FixedField");
            end if;
            ED := EndomorphismData(GeoEndList, L, Kint, Gp, GensHintp, Gphi,
                GeoFactorsQQ, Shorthand : AddTensor := true, AddRing := false,
                AddSatoTate := false, AddDecomposition := false);
            if ED[2] eq ["M_2(RR)"] then
                return "C_{6,1}";
            elif ED[2] eq ["HH"] then
                return "J(C_3)";
            end if;
        elif IsIsomorphic(Hp, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            // In this case it suffices to check whether the polynomial that
            // defines the center of the geometric endomorphism ring in fact
            // has a root in the ground field:
            F := BaseRing(GeoFactorsQQ[1]);
            f := DefiningPolynomial(F);
            if HasRoot(f, K) then
                return "D_2";
            else
                return "J(C_2)";
            end if;
        elif IsIsomorphic(Hp, DirectProduct(CyclicGroup(4), CyclicGroup(2))) then
            return "J(C_4)";
        elif IsIsomorphic(Hp, DirectProduct(CyclicGroup(6), CyclicGroup(2))) then
            return "J(C_6)";
        elif IsIsomorphic(Hp, DihedralGroup(3)) then
            return "D_3";
        elif IsIsomorphic(Hp, DihedralGroup(4)) then
            return "D_4";
        elif IsIsomorphic(Hp, DihedralGroup(6)) then
            return "D_6";
        elif IsIsomorphic(Hp, AlternatingGroup(4)) then
            return "T";
        elif IsIsomorphic(Hp, SymmetricGroup(4)) then
            return "O";
        end if;
    elif FactorsRR eq ["CC", "CC"] then
        Hp := sub<Gp | GensHp>;
        if IsIsomorphic(Hp, CyclicGroup(2)) then
            return "C_2";
        elif IsIsomorphic(Hp, CyclicGroup(3)) then
            return "C_3";
        elif IsIsomorphic(Hp, CyclicGroup(4)) then
            return "C_4";
        elif IsIsomorphic(Hp, CyclicGroup(6)) then
            return "C_6";
        end if;
    elif FactorsRR eq ["M_2(RR)"] then
        return "C_{2,1}";
    elif FactorsRR eq ["HH"] then
        return "J(C_1)";
    elif FactorsRR eq ["M_2(CC)"] then
        return "C_1";
    end if;
end if;

error Error("All cases in SatoTateGroupG2 fell through");

end function;
