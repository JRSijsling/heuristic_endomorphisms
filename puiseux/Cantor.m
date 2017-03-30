/***
 *  Cantor equation functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

forward CantorEquations;

forward CandidateFunctions;
forward FunctionValuesFromApproximations;
forward FunctionsFromApproximations;
forward FunctionsCheck;

forward CantorMorphismFromMatrix;


import "Divisor.m": InitializeCurve;
import "LocalInfo.m": PuiseuxRamificationIndex, InitializeImageBranch;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;


function CantorEquations(f, g);

R := Parent(f);
x,y := Explode(GeneratorsSequence(R));
F := BaseRing(R);
S := PolynomialRing(F, 2*g);
varnames := [ Sprintf("a%o", i) : i in [1..g] ] cat [ Sprintf("b%o", i) : i in [1..g] ];
AssignNames(~S, varnames);
T<t> := PolynomialRing(S);
/* Start with trace and end with norm: */
canpol := t^g + &+[ S.i * t^(g - i) : i in [1..g] ];
substpol := &+[ S.(g + i) * t^(g - i) : i in [1..g] ];
eqpol := Evaluate(f, [t, substpol]) mod canpol;
return Coefficients(eqpol);

end function;


function CandidateFunctions(X, d)

x,y := Explode(GeneratorsSequence(X`R));
f := X`DEs[1];
dens := [ x^i : i in [0..d] ];
nums := [ x^i*y^j : i in [0..d], j in [0..(Degree(f, y) - 1)] | i + j le d ];
return dens, nums;

end function;


function FunctionValuesFromApproximations(X, Qs)

PR := Parent(Qs[1][1]);
R<t> := PolynomialRing(PR);
pol_approx := &*[ t - Q[1] : Q in Qs ];
/* Start with trace and end with norm: */
as_approx := Reverse(Coefficients(pol_approx)[1..X`g]);
v := Matrix([ [ Q[2] : Q in Qs ] ]);
M := Transpose(Matrix([ [ Q[1]^i : i in [0..(X`g - 1)] ] : Q in Qs ]));
w := v*M^(-1);
bs_approx := Reverse(Eltseq(w));
return as_approx cat bs_approx;

end function;


function FunctionsFromApproximations(X, P, Qs, d)

/* Check for debugging: */
//print Valuation(Evaluate(X`DEs[1], P)); print [ Valuation(Evaluate(X`DEs[1], Q)) : Q in Qs ];

I := ideal<X`R | X`DEs[1]>;
dens, nums := CandidateFunctions(X, d);
fs_approx := FunctionValuesFromApproximations(X, Qs);
fs := [ ];
test := true;
for f_approx in fs_approx do
    ev_dens := [ -f_approx * Evaluate(den, P) : den in dens ];
    ev_nums := [ Evaluate(num, P) : num in nums ];
    evs := ev_dens cat ev_nums;
    prec := Floor(Minimum([ AbsolutePrecision(ev) : ev in evs ]));
    M := Matrix([ [ X`F ! Coefficient(ev, i) : i in [0..(prec - 1)] ] : ev in evs ]);
    Ker := Kernel(M);
    if Dimension(Ker) eq 0 then
        test := false;
        fs := [];
        break;
        //Append(~fs, X`K ! 0);
    else
        found := false;
        B := Basis(Ker);
        for b in B do
            v := Eltseq(b);
            /* In general some cancellation takes place here: */
            f := &+[ v[i + #dens]*nums[i] : i in [1..#nums] ] / &+[ v[i]*dens[i] : i in [1..#dens] ];
            //print f;
            //if not (X`R ! Numerator(X`K ! f)) in I then
                Append(~fs, X`K ! f);
                found := true;
                break;
            //end if;
        end for;
        if not found then
            test := false;
            fs := [];
            break;
            //Append(~fs, X`K ! 0);
        end if;
    end if;
end for;
return test, fs;

end function;


function FunctionsCheck(X, fs)

I := ideal<X`R | X`DEs[1]>;
for cantor_eq in X`cantor_eqs do
    if not (X`R ! Numerator(X`K ! Evaluate(cantor_eq, fs))) in I then
        return false;
    end if;
end for;
return true;

end function;


intrinsic CantorMorphismFromMatrix(X::Crv, P0::Pt, M::AlgMatElt : Margin := 2^4, LowerBound := 1) -> Sch
{Given a curve X, a point P0 of X, and a matrix M that gives the tangent
representation of an endomorphism on the standard basis of differentials,
returns a corresponding Cantor morphism (if it exists). The parameter Margin
indicates how many potentially superfluous terms are used in the development of
the branch, and the parameter LowerBound specifies at which degree one starts
to look for a divisor.}

/* We start at a suspected estimate and then increase degree until we find an
 * appropriate divisor: */
e := PuiseuxRamificationIndex(M);
output := InitializeCurve(X, P0);
d := LowerBound;

/* Some global elements needed below: */
g := X`g;
NormM := X`T * M * (X`T)^(-1);

while true do
    vprintf EndoCheck : "Trying degree %o...\n", d;
    dens, nums := CandidateFunctions(X, d);
    n := #dens + #nums + Margin;
    vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

    /* TODO: This does some work many times over.
     * On the other hand, an iterator also has its disadvantages because of superfluous coefficients. */
    vprintf EndoCheck : "Expanding... ";
    P, Qs := ApproximationsFromTangentAction(X, NormM, n*e);
    vprint EndoCheck, 3 : P, Qs;
    vprintf EndoCheck : "done.\n";

    /* Fit a Cantor morphism to it: */
    vprintf EndoCheck : "Solving linear system... ";
    test, fs := FunctionsFromApproximations(X, P, Qs, d);
    vprintf EndoCheck : "done.\n";

    if test then
        as := fs[1..g];
        bs := fs[(g + 1)..(2*g)];
        vprintf EndoCheck : "Checking:\n";
        vprintf EndoCheck : "Step 1... ";
        //test1 := &and[ IsWeaklyZero(Q[1]^g + &+[ Evaluate(as[i], P) * Q[1]^(g - i) : i in [1..g] ]) : Q in Qs ];
        //test2 := &and[ IsWeaklyZero(Q[2]   - &+[ Evaluate(bs[i], P) * Q[1]^(g - i) : i in [1..g] ]) : Q in Qs ];
        //vprintf EndoCheck : "done.\n";
        //if test1 and test2 then
            //vprintf EndoCheck : "Step 2...\n";
            test := FunctionsCheck(X, fs);
            vprintf EndoCheck : "done.\n";
            if test then
                vprintf EndoCheck : "Functions found!\n";
                return fs;
            end if;
        //end if;
    end if;

    /* If that does not work, give up and try one degree higher: */
    d +:= 1;
end while;

end intrinsic;


intrinsic CantorMorphismFromMatrixSplit(X::Crv, P0::Pt, M::AlgMatElt : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), B := 300) -> Sch
{Given a curve X, a point P0 of X, and a matrix M that gives the tangent
representation of an endomorphism on the standard basis of differentials,
returns a corresponding Cantor morphism (if it exists). The parameter Margin
indicates how many potentially superfluous terms are used in the development of
the branch, and the parameter LowerBound specifies at which degree one starts
to look for a divisor.}

/* We start at a suspected estimate and then increase degree until we find an
 * appropriate divisor: */
e := PuiseuxRamificationIndex(M);
output := InitializeCurve(X, P0);
d := LowerBound;
NormM := X`T * M * (X`T)^(-1);
tjs0, f := InitializeImageBranch(M);

/* Some global elements needed below: */
g := X`g;
F := X`F;
rF := X`rF;
OF := X`OF;
BOF := X`BOF;
R := X`R;
K := X`K;
/* TODO: Play with precision here */
P, Qs := ApproximationsFromTangentAction(X, NormM, g);

ps_rts := [ ];
prs := [ ];
I := ideal<X`OF | 1>;
fss_red := [* *];
have_to_check := true;
while true do
    /* Find new prime */
    repeat
        p_rt := RandomSplitPrime(f, B);
        p, rt := Explode(p_rt);
    until not p in [ tup[1] : tup in ps_rts ];
    Append(~ps_rts, p_rt);
    vprintf EndoCheck : "Split prime over %o\n", p;

    /* Add corresponding data: */
    pr := ideal<OF | [ p, rF - rt ]>;
    Append(~prs, pr);
    I *:= pr;
    X_red := ReduceCurveSplit(X, p, rt);
    NormM_red := ReduceMatrixSplit(NormM, p, rt);
    BI := Basis(I);

    /* Uncomment for check on compatibility with reduction */
    //print CantorMorphismFromMatrix(X_red`U, X_red`P0, (X_red`T)^(-1) * M_red * X_red`T);

    while d le UpperBound do
        vprintf EndoCheck : "Trying degree %o...\n", d;
        dens_red, nums_red := CandidateFunctions(X_red, d);
        n := #dens_red + #nums_red + Margin;
        vprintf EndoCheck : "Number of digits in expansion: %o.\n", n*e;

        /* Take non-zero image branch: */
        /* TODO: This does some work many times over, but only the first time */
        vprintf EndoCheck, 2 : "Expanding... ";
        P_red, Qs_red := ApproximationsFromTangentAction(X_red, NormM_red, n*e);
        vprint EndoCheck, 3 : P_red, Qs_red;
        vprintf EndoCheck, 2 : "done.\n";

        /* Fit a Cantor morphism to it: */
        vprintf EndoCheck, 2 : "Solving linear system... ";
        test_red, fs_red := FunctionsFromApproximations(X_red, P_red, Qs_red, d);
        vprintf EndoCheck, 2 : "done.\n";

        if test_red then
            vprintf EndoCheck, 2 : "Checking:\n";
            vprintf EndoCheck, 2 : "Step 1: ";
            if not have_to_check then
                vprintf EndoCheck, 2 : "done.\n";
                vprintf EndoCheck, 2 : "Functions found!\n";
                break;
            end if;
            test := FunctionsCheck(X_red, fs_red);
            vprintf EndoCheck, 2 : "done.\n";
            if test then
                vprintf EndoCheck, 2 : "Functions found!\n";
                have_to_check := false;
                break;
            end if;
        end if;

        /* If that does not work, give up and try one degree higher. */
        d +:= 1;
    end while;
    Append(~fss_red, fs_red);

    vprintf EndoCheck : "Fractional CRT... ";
    fs := [ ];
    for i:=1 to #fss_red[1] do
        num := R ! 0;
        for mon in Monomials(Numerator(fss_red[1][i])) do
            exp := Exponents(mon);
            rs := [* MonomialCoefficient(Numerator(fss_red[j][i]), exp) : j in [1..#fss_red] *];
            num +:= FractionalCRTSplit(rs, prs, OF, I, BOF, BI, F) * Monomial(R, exp);
        end for;
        den := R ! 0;
        for mon in Monomials(Denominator(fss_red[1][i])) do
            exp := Exponents(mon);
            rs := [* MonomialCoefficient(Denominator(fss_red[j][i]), exp) : j in [1..#fss_red] *];
            den +:= FractionalCRTSplit(rs, prs, OF, I, BOF, BI, F) * Monomial(R, exp);
        end for;
        Append(~fs, K ! (num / den));
    end for;
    vprintf EndoCheck : "done.\n";

    if #fs ne 2*g then
        print "Algorithm failed to find a Cantor morphism";
        return [ ];
    end if;
    vprintf EndoCheck : "Checking:\n";
    vprintf EndoCheck : "Step 1... ";
    as := fs[1..g];
    bs := fs[(g + 1)..(2*g)];
    test1 := true;
    for Q in Qs do
        if not IsWeaklyZero(Q[1]^g + &+[ Evaluate(as[i], P) * Q[1]^(g - i) : i in [1..g] ]) then
            test1 := false;
            break;
        end if;
        if not IsWeaklyZero(Q[2]   - &+[ Evaluate(bs[i], P) * Q[1]^(g - i) : i in [1..g] ]) then
            test1 := false;
            break;
        end if;
    end for;
    vprintf EndoCheck : "done.\n";

    if test1 then
        vprintf EndoCheck : "Step 2...\n";
        vprintf EndoCheck : "Candidate functions:\n";
        vprint EndoCheck : fs;
        test2 := FunctionsCheck(X, fs);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Functions found!\n";
            return fs;
        end if;
    end if;
end while;

end intrinsic;

