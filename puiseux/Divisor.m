/***
 *  Divisor functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


declare attributes Crv : is_hyperelliptic, is_planar, is_smooth, is_plane_quartic;
declare attributes Crv : unif, unif_index;
declare attributes Crv : g, U, P0, A, DEs;
declare attributes Crv : patch_index, R, x, y, K;
declare attributes Crv : F, rF, OF, BOF;
declare attributes Crv : OurB, NormB, T;
declare attributes Crv : initialized;
declare attributes Crv : cantor_equations;

declare verbose EndoCheck, 3;


forward InitializeCurve;
forward AlgebraicUniformizer;
forward OurBasisOfDifferentials;
forward ChangePatchTangentAction;
forward NormalizedBasisOfDifferentials;

forward VariableOrder;
forward ExtractPoint;
forward ExtractPoints;
forward ExtractHomomorphisms;

forward CandidateDivisors;
forward IrreducibleComponentsFromBranches;
forward CheckEquations;
forward CheckIrreducibleComponent;

forward DivisorFromMatrixSplit;
forward DivisorFromMatrixSplitStepModP;


import "LocalInfo.m": DevelopPoint, InitializeImageBranch;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Cantor.m": CantorEquations;


function InitializeCurve(X, P0)
/* This is written in a silly way; it is essentially a procedure */

if not assigned X`initialized then
    X`initialized := false;
end if;
if X`initialized then
    return 0;
end if;

X`is_hyperelliptic := IsHyperelliptic(X); X`is_planar := IsPlaneCurve(X); X`is_smooth := IsNonSingular(X);
X`g := Genus(X); X`is_plane_quartic := (X`is_planar) and (X`is_smooth) and (X`g eq 3);

/*Find affine patch and coordinates */
if IsAffine(X) then
    X`U := X; X`P0 := P0; X`patch_index := 1;
else
    X`U, X`P0, X`patch_index := AffinePatch(X, P0);
end if;
X`A := Ambient(X`U); X`R := CoordinateRing(X`A);

/* If the equation comes in the form y^2 + h y = f, or a transformed version of
 * this, then we name x and y in that way; otherwise we just use the ambient.
 */
if (X`is_hyperelliptic or X`g eq 1) and (X`patch_index eq 3) then
    X`x := X`R.2; X`y := X`R.1;
else
    X`x := X`R.1; X`y := X`R.2;
end if;

/* Construct equation order */
X`F := BaseRing(X`R);
if Type(X`F) eq FldRat then
    X`rF := 1;
    X`OF := Integers();
else
    X`rF := Denominator(X`F.1) * X`F.1;
    X`OF := Order([ X`rF^i : i in [0..Degree(X`F) - 1] ]);
end if;

X`BOF := Basis(X`OF); X`K := FieldOfFractions(X`R);
X`DEs := DefiningEquations(X`U);
X`unif, X`unif_index := AlgebraicUniformizer(X);

X`OurB := OurBasisOfDifferentials(X);
X`NormB, X`T := NormalizedBasisOfDifferentials(X);

// TODO: Perhaps this steps costs a bit. Not removed for now.
if X`is_planar then
    X`cantor_equations := CantorEquations(X);
end if;
X`initialized := true;
return 0;

end function;


function AlgebraicUniformizer(X)
/*
 * Input:   A plane curve X.
 * Output:  A uniformizing element at P0
 *          and the corresponding index.
 */

fX := X`DEs[1]; x := X`x; y := X`y;
/* We use x when possible, and then return it along with its index */
if Evaluate(Derivative(fX, y), X`P0) ne 0 then
    if X`R.1 eq x then
        return x, 1;
    else
        return x, 2;
    end if;
end if;

/* Otherwise we use y */
if X`R.1 eq y then
    return y, 1;
else
    return y, 2;
end if;

end function;


function OurBasisOfDifferentials(X)
/*
 * Input:   A curve X.
 * Output:  A basis of global differentials on X, represented by elements of
 *          the rational function field by using our choice of uniformizer
 */

g := X`g;
R := X`R; x := X`x; y := X`y;
f := X`DEs[1];

if g eq 0 then
    return [ ];

elif X`is_hyperelliptic or (g eq 1) then
    /* (Hyper)elliptic case: we use x^i dx / y */
    s := MonomialCoefficient(f, y^2);
    if g eq 1 then
        if X`unif eq x then
            return [  2*s / Derivative(f, y) ];
        else
            return [ -2*s / Derivative(f, x) ];
        end if;
    end if;
    /* Otherwise the uniformizer will always be x, since we have to use a
     * non-Weierstrass point */
    return [ 2*s*x^(i-1) / Derivative(f, y) : i in [1..g] ];

elif X`is_plane_quartic then
    /* Plane quartic case: we use ({x,y,1} / (dF / dy)) dx */
    if X`unif_index eq 1 then
        return [ X`K ! ( n / Derivative(f, 2)) : n in [R.1, R.2, 1] ];
    else
        return [ X`K ! (-n / Derivative(f, 1)) : n in [R.1, R.2, 1] ];
    end if;

else
    error "OurBasisOfDifferentials not implemented yet for this curve";
end if;

end function;


function ChangePatchTangentAction(X, Y, M)
/*
 * Input:  Two curves X, Y and a representation M on the standard basis of
 *         differentials.
 * Output: Matrix for standard differentials on the patches of X and Y used.
 */

/* M acts on the right, so to precompose with the operation on X we multiply on
 * the left; we modify the rows. */
if X`g eq 1 or X`is_hyperelliptic then
    /* In this case only patch index 3 is possible */
    if X`patch_index eq 3 then
        M := Matrix(X`F, [ Reverse([ -c : c in Eltseq(row)]) : row in Rows(M) ]);
    end if;

elif X`is_plane_quartic then
    if X`patch_index eq 2 then
        M := Matrix(X`F, [ [ -row[1], -row[3], -row[2] ] : row in Rows(M) ]);
    elif X`patch_index eq 3 then
        M := Matrix(X`F, [ [ row[2], row[3], row[1] ] : row in Rows(M) ]);
    end if;
end if;

/* For y we postcompose, hence we have to modify columns; we therefore take a
 * transpose and go back */
M := Transpose(M);
if Y`g eq 1 or Y`is_hyperelliptic then
    /* In this case only patch index 3 is possible */
    if Y`patch_index eq 3 then
        M := Matrix(Y`F, [ Reverse([ -c : c in Eltseq(row)]) : row in Rows(M) ]);
    end if;

elif Y`is_plane_quartic then
    if Y`patch_index eq 2 then
        M := Matrix(Y`F, [ [ -row[1], -row[3], -row[2] ] : row in Rows(M) ]);
    elif Y`patch_index eq 3 then
        M := Matrix(Y`F, [ [ row[2], row[3], row[1] ] : row in Rows(M) ]);
    end if;
end if;

M := Transpose(M);
return M;

end function;


function NormalizedBasisOfDifferentials(X)
/*
 * Input:   A curve X.
 * Output:  A differential basis Bnorm that is normalized with respect to the uniformizing parameter,
 *          and a matrix T such that multiplication by T on the left sends B to Bnorm.
 */

P := DevelopPoint(X, X`P0, X`g);
BP := [ Evaluate(b, P) : b in X`OurB ];
T := Matrix([ [ Coefficient(BP[i], j - 1) : j in [1..X`g] ] : i in [1..X`g] ])^(-1);
NormB := [ &+[ T[i,j] * X`OurB[j] : j in [1..X`g] ] : i in [1..X`g] ];
return NormB, T;

end function;


function VariableOrder()
/*
 * The order in which x(P), y(P), x(Q), y(Q) and hence x1, y1, x2, y2 are used
 * in the product space. Note that because of the lexicographical ordering
 * variables that occur later are eliminated for first.
 */

/* x(P) to 4th comp, y(P) to 2nd comp, etc */
return [4, 2, 3, 1];

end function;


function ExtractPoint(X, P)
/* Finds x(P) and y(P); indices depend on the choice of patch */

R := X`R; x := X`x; y := X`y;
if x eq R.1 then
    return [ P[1], P[2] ];
else
    return [ P[2], P[1] ];
end if;

end function;


function ExtractPoints(X, Y, P, Q)
/* Because x and y may change depending on the choice of patch */

seq := ExtractPoint(X, P) cat ExtractPoint(Y, Q);
varord := VariableOrder();
return [ seq[varord[i]] : i in [1..#varord] ];

end function;


function ExtractHomomorphisms(X, Y)

RX := X`R; RY := Y`R;
varord := VariableOrder();
Rprod := PolynomialRing(X`F, 4, "lex");
seqX := [ Rprod.varord[i] : i in [1..2] ];
seqY := [ Rprod.varord[i] : i in [3..4] ];
if X`x ne RX.1 then
    Reverse(~seqX);
end if;
if Y`x ne RY.1 then
    Reverse(~seqY);
end if;
hX := hom<RX -> Rprod | seqX >; hY := hom<RY -> Rprod | seqY >;
return [ hX, hY ];

end function;


function CandidateDivisors(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

gX := X`g; fX := X`DEs[1]; RX := X`R; xX := X`x; yX := X`y;
gY := Y`g; fY := Y`DEs[1]; RY := Y`R; xY := Y`x; yY := Y`y;

if X`is_hyperelliptic then
    divsX := [ xX^i : i in [0..(d div 2)] ] cat [ xX^i*yX : i in [0..((d - gX - 1) div 2)] ];
    divsY := [ xY^i : i in [0..(d div 2)] ] cat [ xY^i*yY : i in [0..((d - gX - 1) div 2)] ];
    divsY := [ xY^i : i in [0..gY] ] cat [ yY ];
    Reverse(~divsX); Reverse(~divsY);
elif X`is_planar then
    divsX := [ xX^i*yX^j : i in [0..d] , j in [0..(Degree(fX, yX) - 1)] | i + j le d  ];
    divsY := [ xY^i*yY^j : i in [0..gY], j in [0..(Degree(fY, yY) - 1)] | i + j le gY ];
    //divsY := [ xY^i : i in [0..gY] ] cat [ yY ];
    Reverse(~divsX); Reverse(~divsY);
end if;

varord := VariableOrder();
hs := ExtractHomomorphisms(X, Y);
CP := [ [* divX, divY *] : divX in divsX, divY in divsY ];
divs := [ &*[ hs[i](tup[i]) : i in [1..2] ] : tup in CP ];
divs := Reverse(Sort(divs));
return divs;

end function;


function IrreducibleComponentsFromBranches(X, Y, fs, P, Qs : DivPP1 := false)
/*
 * Input:   Two curves X and Y,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components that fit the given data.
 */

/* Recovering a linear system */
e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
prec := Precision(Parent(P[1]));
M := [ ];
for f in fs do
    r := [ ];
    for Q in Qs do
        seq := ExtractPoints(X, Y, P, Q);
        ev := Evaluate(f, seq);
        r cat:= [ Coefficient(ev, i/e) : i in [0..prec - X`g] ];
    end for;
    Append(~M, r);
end for;
B := Basis(Kernel(Matrix(M)));
/* Coerce back to ground field (possible because of echelon form) */
B := [ [ X`F ! c : c in Eltseq(b) ] : b in B ];

/* Corresponding equations */
hX, hY := Explode(ExtractHomomorphisms(X, Y));
Rprod := Codomain(hX);
eqs := [ Rprod ! (&+[ b[i] * fs[i] : i in [1..#fs] ]) : b in B ];
eqs := eqs cat [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ];

if DivPP1 then
    vprintf EndoCheck, 2 : "Calculating final element in Groebner basis... ";
    Rprod := Codomain(hX);
    GB := GroebnerBasis(ideal< Rprod | eqs >);
    vprint EndoCheck, 3 : GB;
    vprintf EndoCheck, 2 : "done.\n";
    Append(~eqs, GB[#GB]);
end if;

/* Corresponding scheme */
A := AffineSpace(Rprod);
S := Scheme(A, eqs);
return [ S ];

/* TODO: These steps may be a time sink and should be redundant, so we avoid
 *       them. They get eliminated as the degree increases anyway. */
return [ ReducedSubscheme(I) : I in IrreducibleComponents(S) ];

end function;


function CheckEquations(X, Y, P, Qs, DEs)

for DE in DEs do
    for Q in Qs do
        seq := ExtractPoints(X, Y, P, Q);
        if not IsWeaklyZero(Evaluate(DE, seq)) then
            return false;
        end if;
    end for;
end for;
return true;

end function;


function CheckIrreducibleComponent(X, Y, I)
/*
 * Input:   An irreducible scheme I in X x Y.
 * Output:  Whether or not I intersects P0 x X with the correct multiplicity at
 *          P0 and nowhere else.
 */

A4 := Ambient(I); R4 := CoordinateRing(A4);
R2 := PolynomialRing(X`F, 2); A2 := AffineSpace(R2);
/* Note that here we do not care about the order of x and y */
seq := [ X`P0[1], X`P0[2], R2.1, R2.2 ];
varord := VariableOrder();
seq := [ seq[varord[i]] : i in [1..#varord] ];
h := hom< R4 -> R2 | seq >;
eqs2 := [ h(eq4) : eq4 in DefiningEquations(I) ];
S := Scheme(A2, eqs2);
if Dimension(S) eq 0 then
    if Degree(ReducedSubscheme(S)) eq 1 then
        if Degree(S) eq Y`g then
            //if Dimension(I) eq 1 then
                return true;
            //end if;
        end if;
    end if;
end if;
return false;

end function;


intrinsic DivisorFromMatrix(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), DivPP1 := false) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);
NormM := ChangePatchTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);

d := LowerBound;
while true do
    vprintf EndoCheck : "Trying degree %o...\n", d;
    fs := CandidateDivisors(X, Y, d);
    n := #fs + Margin;
    vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

    /* Take non-zero image branch */
    vprintf EndoCheck : "Expanding... ";
    P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n);
    vprint EndoCheck, 3 : P, Qs;
    vprintf EndoCheck : "done.\n";

    /* Fit a divisor to it */
    vprintf EndoCheck : "Solving linear system... ";
    ICs := IrreducibleComponentsFromBranches(X, Y, fs, P, Qs : DivPP1 := DivPP1);
    vprintf EndoCheck : "done.\n";

    for S in ICs do
        DEs := DefiningEquations(S);
        vprintf EndoCheck : "Checking:\n";
        vprintf EndoCheck : "Step 1... ";
        test1 := CheckEquations(X, Y, P, Qs, DEs);
        vprintf EndoCheck : "done.\n";
        if test1 then
            vprintf EndoCheck : "Step 2... ";
            test2 := CheckIrreducibleComponent(X, Y, S);
            vprintf EndoCheck : "done.\n";
            if test2 then
                vprintf EndoCheck : "Divisor found!\n";
                return true, S;
            end if;
        end if;
    end for;

    /* If that does not work, give up and try one degree higher */
    d +:= 1;
    if d gt UpperBound then
        return false, "";
    end if;
end while;

end intrinsic;


intrinsic DivisorFromMatrixSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), DivPP1 := false, B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor */
output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);
NormM := ChangePatchTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);
tjs0, f := InitializeImageBranch(NormM);
F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF;
Rprod := PolynomialRing(X`F, 4, "lex");
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, X`g);

ps_rts := [ ]; prs := [ ]; DEss_red := [* *];
I := ideal<X`OF | 1>;
have_to_check := true;

d := LowerBound;
while true do
    /* Find new prime */
    repeat
        p_rt := RandomSplitPrime(f, B);
        p, rt := Explode(p_rt);
    until not p in [ tup[1] : tup in ps_rts ];
    Append(~ps_rts, p_rt);
    vprintf EndoCheck : "Split prime over %o\n", p;

    /* Add corresponding data */
    pr := ideal<X`OF | [ p, rF - rt ]>;
    Append(~prs, pr); I *:= pr;
    X_red := ReduceCurveSplit(X, p, rt); Y_red := ReduceCurveSplit(Y, p, rt);
    NormM_red := ReduceMatrixSplit(NormM, p, rt);
    BI := Basis(I);
    while true do
        found, S_red := DivisorFromMatrixSplitStepModP(X_red, Y_red, NormM_red, d : Margin := Margin, DivPP1 := DivPP1, have_to_check := have_to_check);
        /* If that does not work, give up and try one degree higher. Note that
         * d is initialized in the outer loop, so that we keep the degree that
         * works. */
        if found then
            break;
        end if;
        d +:= 1;
        if d gt UpperBound then
            return false, "";
        end if;
    end while;
    have_to_check := false;
    Append(~DEss_red, DefiningEquations(S_red));

    vprintf EndoCheck : "Fractional CRT... ";
    DEs := [ ];
    for i:=1 to #DEss_red[1] do
        DE := Rprod ! 0;
        for mon in Monomials(DEss_red[1][i]) do
            exp := Exponents(mon);
            rs := [* *];
            for j:=1 to #DEss_red do
                Rprod_red := Parent(DEss_red[j][1]);
                Append(~rs, MonomialCoefficient(DEss_red[j][i], Monomial(Rprod_red, exp)));
            end for;
            DE +:= FractionalCRTSplit(rs, prs, OF, I, BOF, BI, F) * Monomial(Rprod, exp);
        end for;
        Append(~DEs, DE);
    end for;
    vprintf EndoCheck : "done.\n";

    vprintf EndoCheck : "Checking:\n";
    vprintf EndoCheck : "Step 1... ";
    /* Note that P and Qs are calculated at the beginning of this function */
    test1 := CheckEquations(X, Y, P, Qs, DEs);
    vprintf EndoCheck : "done.\n";

    if test1 then
        S := Scheme(AffineSpace(Rprod), DEs);
        vprintf EndoCheck : "Step 2...\n";
        test2 := CheckIrreducibleComponent(X, Y, S);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Divisor found!\n";
            return true, S;
        end if;
    end if;
end while;

end intrinsic;


function DivisorFromMatrixSplitStepModP(X_red, Y_red, NormM_red, d : Margin := 2^4, DivPP1 := false, have_to_check := true)
/* Step mod p of the above */

vprintf EndoCheck : "Trying degree %o...\n", d;
fs_red := CandidateDivisors(X_red, Y_red, d);
n := #fs_red + Margin;
vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

/* Take non-zero image branch */
vprintf EndoCheck, 2 : "Expanding... ";
P_red, Qs_red := ApproximationsFromTangentAction(X_red, Y_red, NormM_red, n);
vprint EndoCheck, 3 : P_red, Qs_red;
vprintf EndoCheck, 2 : "done.\n";

/* Fit a divisor to it */
vprintf EndoCheck, 2 : "Solving linear system... ";
ICs_red := IrreducibleComponentsFromBranches(X_red, Y_red, fs_red, P_red, Qs_red : DivPP1 := DivPP1);
vprintf EndoCheck, 2 : "done.\n";

for S_red in ICs_red do
    vprintf EndoCheck, 2 : "Checking:\n";
    vprintf EndoCheck, 2 : "Step 1... ";

    if not have_to_check then
        vprintf EndoCheck, 2 : "done.\n";
        vprintf EndoCheck, 2 : "Divisor found!\n";
        return true, S_red;
    end if;

    test := CheckIrreducibleComponent(X_red, Y_red, S_red);
    vprintf EndoCheck, 2 : "done.\n";
    if test then
        vprintf EndoCheck, 2 : "Divisor found!\n";
        return true, S_red;
    end if;
end for;
return false, 0;

end function;
