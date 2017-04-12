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


declare attributes Crv : is_hyp, is_planar, is_smooth, is_plane_quartic;
declare attributes Crv : unif, unif_index;
declare attributes Crv : g, U, P0, A, R, F, rF, OF, BOF, K, DEs, OurB, NormB, T;
declare attributes Crv : initialized;
declare attributes Crv : cantor_eqs;

declare verbose EndoCheck, 3;


forward InitializeCurve;
forward AlgebraicUniformizer;
forward OurBasisOfDifferentials;
forward NormalizedBasisOfDifferentials;

forward CandidateDivisors;
forward IrreducibleComponentsFromBranches;
forward IrreducibleComponentCheck;

forward DivisorFromMatrix;
forward DivisorFromMatrixSplit;


import "LocalInfo.m": DevelopPoint, InitializeImageBranch;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Cantor.m": CantorEquations;


function InitializeCurve(X, P0)

if not assigned X`initialized then
    X`initialized := false;
end if;
if X`initialized then
    return 0;
end if;
X`is_hyp := IsHyperelliptic(X); X`is_planar := IsPlaneCurve(X); X`is_smooth := IsNonSingular(X);
X`g := Genus(X); X`is_plane_quartic := (X`is_planar) and (X`is_smooth) and (X`g eq 3);
if IsAffine(X) then
    X`U := X; X`P0 := P0;
else
    X`U, X`P0 := AffinePatch(X, P0);
end if;
X`A := Ambient(X`U); X`R := CoordinateRing(X`A); X`F := BaseRing(X`R);
if Type(X`F) eq FldRat then
    X`rF := 1; X`OF := Integers();
else
    X`rF := Denominator(X`F.1) * X`F.1;
    X`OF := Order([ X`rF^i : i in [0..Degree(X`F) - 1] ]);
end if;
X`BOF := Basis(X`OF); X`K := FieldOfFractions(X`R);
X`DEs := DefiningEquations(X`U);
X`unif, X`unif_index := AlgebraicUniformizer(X);
X`OurB := OurBasisOfDifferentials(X);
X`NormB, X`T := NormalizedBasisOfDifferentials(X);
if X`is_planar then
    X`cantor_eqs := CantorEquations(X`DEs[1], X`g);
end if;
X`initialized := true;
return 0;

end function;


function AlgebraicUniformizer(X)
/*
 * Input:   A curve X.
 * Output:  A uniformizing element at P0
 *          and the corresponding index.
 */

Gens := GeneratorsSequence(X`R);
M := Matrix([ [ Evaluate(Derivative(DE, gen), X`P0) : gen in Gens ] : DE in X`DEs ]);
/* Default is the first coordinate: */
i0 := 1;
for i in [1..#Gens] do
    if &and[ M[j, i] eq 0 : j in [1..#Rows(M)] ] then
        i0 := i;
    end if;
end for;
return Gens[i0], i0;

end function;


function OurBasisOfDifferentials(X);
/*
 * Input:   A curve X.
 * Output:  A basis of global differentials on X, represented by elements of
 *          the ambient.
 */

g := X`g;
R := X`R;
x := R.1; y := R.2;
if X`is_hyp then
    f := X`DEs[1];
    c2 := MonomialCoefficient(f, y^2);
    c1 := &+[ MonomialCoefficient(f, x^i*y) * x^i : i in [0..g] ];
    return [ x^(i-1) / (y + c1/(2*c2)) : i in [1..g] ];
elif X`is_plane_quartic then
    f := X`DEs[1];
    if X`unif_index eq 1 then
        return [ X`K ! (n / Derivative(f, 2)) : n in [x,y,1] ];
    else
        return [ X`K ! (n / Derivative(f, 1)) : n in [x,y,1] ];
    end if;
else
    B := BasisOfDifferentialsFirstKind(X`U);
    du := Differential(AlgebraicUniformizer(X));
    return [ X`K ! (b / du) : b in B ];
end if;

end function;


function NormalizedBasisOfDifferentials(X)
/*
 * Input:   A curve X.
 * Output:  A differential basis Bnorm that is normalized with respect to the uniformizing parameter,
 *          and a matrix T such that multiplication by T sends B to Bnorm.
 */

g := X`g;
P := DevelopPoint(X, X`P0, g);
BP := [ Evaluate(b, P) : b in X`OurB ];
T := Matrix([ [ Coefficient(BP[i], j - 1) : j in [1..g] ] : i in [1..g] ])^(-1);
NormB := [ &+[ T[i,j] * X`OurB[j] : j in [1..g] ] : i in [1..g] ];
return NormB, T;

end function;


function CandidateDivisors(X, d)
/*
 * Input:   A curve X
 *          and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

g := X`g; R := X`R; F := X`F;
dim := Rank(R);

/* TODO: Asymmetry */
f := X`DEs[1];
x,y := Explode(GeneratorsSequence(R));
Rprod := PolynomialRing(F, 2 * dim);
if X`is_hyp then
    Xdivs1 := [ x^i : i in [0..(d div 2)] ] cat [ x^i*y : i in [0..((d - g - 1) div 2)] ];
elif X`is_planar then
    f := DefiningEquations(X`U)[1];
    Xdivs1 := [ x^i*y^j : i in [0..d], j in [0..(Degree(f) - 1)] | i + j le d ];
end if;
Xdivs2 := [ x^i*y^j : i in [0..g], j in [0..(Degree(f, 2) - 1)] ];
Xdivs2 := [ x^i : i in [0..g] ] cat [ y ];

hs := [ hom<R -> Rprod | [ Rprod.j : j in [ ((i-1)*dim + 1)..i*dim ] ]> : i in [1..2] ];
CP := CartesianProduct([Xdivs1, Xdivs2]);
return [ &*[ hs[i](tup[i]) : i in [1..2] ] : tup in CP ];

end function;


function IrreducibleComponentsFromBranches(X, fs, P, Qs)
/*
 * Input:   A curve X,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components corresponding that fit the given data.
 */

/* Check for debugging: */
//print Valuation(Evaluate(X`DEs[1], P)); print [ Valuation(Evaluate(X`DEs[1], Q)) : Q in Qs ];

/* Recovering a linear system: */
e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
prec := Precision(Parent(P[1]));
M := [ ];
for f in fs do
    r := [ ];
    for Q in Qs do
        ev := Evaluate(f, P cat Q);
        r cat:= [ Coefficient(ev, i/e) : i in [0..prec - X`g] ];
    end for;
    Append(~M, r);
end for;
M := Matrix(M);
B := Basis(Kernel(M));

/* Coerce back to ground field (possible because of echelon form): */
F := BaseRing(X`U);
B := [ [ F ! c : c in Eltseq(b) ] : b in B ];

/* Corresponding equations: */
DEs := X`DEs; R := X`R; Rprod := Parent(fs[1]); d := Rank(R); g := X`g;
hs := [ hom<R -> Rprod | [ Rprod.j : j in [ ((i-1)*d + 1)..i*d ] ]> : i in [1..2] ];
eqs := [ &+[ b[i] * fs[i] : i in [1..#fs] ] : b in B ];
eqs := eqs cat [ h(DE) : h in hs, DE in DEs ];

/* Corresponding scheme: */
A := AffineSpace(Rprod);
S := Scheme(A, eqs);

return [ S ];

/* TODO: These steps may be a time sink and should be redundant, so we avoid
 *       them. They get eliminated as the degree increases anyway. */
Is := IrreducibleComponents(S);
return [ ReducedSubscheme(I) : I in Is ];

end function;


function IrreducibleComponentCheck(X, I)
/*
 * Input:   An irreducible scheme I in X x X.
 * Output:  Whether or not I intersects P0 x X with the correct multiplicity at
 *          P0 and nowhere else.
 */

A4 := Ambient(I);
R4 := CoordinateRing(A4);
R2 := PolynomialRing(BaseRing(R4), 2);
A2 := AffineSpace(R2);
h := hom< R4 -> R2 | [ X`P0[i] : i in [1..2] ] cat [ R2.i : i in [1..2] ] >;
eqs2 := [ h(eq4) : eq4 in DefiningEquations(I) ];
S := Scheme(A2, eqs2);
if Dimension(S) eq 0 then
    if Degree(ReducedSubscheme(S)) eq 1 then
        if Degree(S) eq X`g then
            /* TODO: This is potentially slightly unsafe but delivers a big speedup */
            //if Dimension(I) eq 1 then
                return true;
            //end if;
        end if;
    end if;
end if;
return false;

end function;


intrinsic DivisorFromMatrix(X::Crv, P0::Pt, M::AlgMatElt : Margin := 2^4, LowerBound := 1) -> Sch
{Given a curve X, a point P0 of X, and a matrix M that gives the tangent
representation of an endomorphism on the standard basis of differentials,
returns a corresponding divisor (if it exists). The parameter Margin indicates
how many potentially superfluous terms are used in the development of the
branch, and the parameter LowerBound specifies at which degree one starts to
look for a divisor.}

/* We start at a suspected estimate and then increase degree until we find an
 * appropriate divisor: */
output := InitializeCurve(X, P0);
d := LowerBound;
NormM := X`T * M * (X`T)^(-1);
while true do
    vprintf EndoCheck : "Trying degree %o...\n", d;
    fs := CandidateDivisors(X, d);
    n := #fs + Margin;
    vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

    /* Take non-zero image branch: */
    vprintf EndoCheck : "Expanding... ";
    P, Qs := ApproximationsFromTangentAction(X, NormM, n);
    vprint EndoCheck, 3 : P, Qs;
    vprintf EndoCheck : "done.\n";

    /* Fit a divisor to it: */
    vprintf EndoCheck : "Solving linear system... ";
    ICs := IrreducibleComponentsFromBranches(X, fs, P, Qs);
    vprintf EndoCheck : "done.\n";

    for S in ICs do
        DEs := DefiningEquations(S);
        vprintf EndoCheck : "Checking:\n";
        vprintf EndoCheck : "Step 1... ";
        //test1 := &and[ &and[ IsWeaklyZero(Evaluate(DE, P cat Q)) : Q in Qs ] : DE in DEs ];
        //vprintf EndoCheck : "done.\n";
        //if test1 then
            //vprintf EndoCheck : "Step 2... ";
            test2 := IrreducibleComponentCheck(X, S);
            vprintf EndoCheck : "done.\n";
            if test2 then
                vprintf EndoCheck : "Divisor found!\n";
                return S;
            end if;
        //end if;
    end for;

    /* If that does not work, give up and try one degree higher: */
    d +:= 1;
end while;

end intrinsic;


intrinsic DivisorFromMatrixSplit(X::Crv, P0::Pt, M::AlgMatElt : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), B := 300) -> Sch
{Given a curve X, a point P0 of X, and a matrix M that gives the tangent
representation of an endomorphism on the standard basis of differentials,
returns a corresponding divisor (if it exists). The parameter Margin indicates
how many potentially superfluous terms are used in the development of the
branch, and the parameter LowerBound specifies at which degree one starts to
look for a divisor.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor: */
output := InitializeCurve(X, P0);
d := LowerBound;
M := X`T * M * (X`T)^(-1);
tjs0, f := InitializeImageBranch(M);

/* Some global elements needed below: */
F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF;
/* TODO: Play with precision here */
P, Qs := ApproximationsFromTangentAction(X, M, X`g);
Rprod := PolynomialRing(X`F, 2 * Rank(X`R));

ps_rts := [ ]; prs := [ ]; I := ideal<X`OF | 1>; DEss_red := [* *]; have_to_check := true;
while true do
    /* Find new prime */
    repeat
        p_rt := RandomSplitPrime(f, B);
        p, rt := Explode(p_rt);
    until not p in [ tup[1] : tup in ps_rts ];
    Append(~ps_rts, p_rt);
    vprintf EndoCheck : "Split prime over %o\n", p;

    /* Add corresponding data: */
    pr := ideal<X`OF | [ p, rF - rt ]>; Append(~prs, pr); I *:= pr;
    X_red := ReduceCurveSplit(X, p, rt); M_red := ReduceMatrixSplit(M, p, rt);
    BI := Basis(I);

    /* Uncomment for check on compatibility with reduction */
    //print DivisorFromMatrix(X_red`U, X_red`P0, (X_red`T)^(-1) * M_red * X_red`T);

    done := false;
    while d le UpperBound do
        vprintf EndoCheck : "Trying degree %o...\n", d;
        fs_red := CandidateDivisors(X_red, d);
        n := #fs_red + Margin;
        vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

        /* Take non-zero image branch: */
        vprintf EndoCheck, 2 : "Expanding... ";
        P_red, Qs_red := ApproximationsFromTangentAction(X_red, M_red, n);
        vprint EndoCheck, 3 : P_red, Qs_red;
        vprintf EndoCheck, 2 : "done.\n";

        /* Fit a divisor to it: */
        vprintf EndoCheck, 2 : "Solving linear system... ";
        ICs_red := IrreducibleComponentsFromBranches(X_red, fs_red, P_red, Qs_red);
        vprintf EndoCheck, 2 : "done.\n";

        for S_red_it in ICs_red do
            vprintf EndoCheck, 2 : "Checking:\n";
            vprintf EndoCheck, 2 : "Step 1... ";
            if not have_to_check then
                vprintf EndoCheck, 2 : "done.\n";
                vprintf EndoCheck, 2 : "Divisor found!\n";
                done := true;
                S_red := S_red_it;
                break;
            end if;
            test := IrreducibleComponentCheck(X_red, S_red_it);
            vprintf EndoCheck, 2 : "done.\n";
            if test then
                vprintf EndoCheck, 2 : "Divisor found!\n";
                done := true;
                S_red := S_red_it;
                have_to_check := false;
                break;
            end if;
        end for;

        if done then
            break;
        end if;

        /* If that does not work, give up and try one degree higher.
         * Note that d is initialized in the outer loop,
         * so that we keep the degree that works. */
        d +:= 1;
    end while;
    if not assigned S_red then
        print "Algorithm failed to find a divisor";
        return 0;
    end if;
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
    test1 := true;
    for DE in DEs do
        test1_int := true;
        for Q in Qs do
            if not IsWeaklyZero(Evaluate(DE, P cat Q)) then
                test1_int := false;
                break;
            end if;
        end for;
        if not test1_int then
            test1 := false;
            break;
        end if;
    end for;
    vprintf EndoCheck : "done.\n";

    if test1 then
        S := Scheme(AffineSpace(Rprod), DEs);
        vprintf EndoCheck : "Step 2...\n";
        vprintf EndoCheck : "Candidate divisor:\n";
        vprint EndoCheck : S;
        test2 := IrreducibleComponentCheck(X, S);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Divisor found!\n";
            return S;
        end if;
    end if;
end while;

end intrinsic;
