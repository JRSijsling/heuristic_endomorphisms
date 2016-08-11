function MyUniformizer(X, P0, is_planar);
/*
 * Input:   A curve X,
 *          a point P0 on X,
 *          and whether or not X is planar.
 * Output:  A uniformizing element at P0.
 */

if is_planar then
    K := FunctionField(X);
    f := DefiningEquations(X)[1];
    R<x, y> := Parent(f);
    if Evaluate(Derivative(f, y), P0) ne 0 then
        return K ! x;
    else
        return K ! y;
    end if;
else
    return UniformizingParameter(P0);
end if;

end function;


function DevelopInUnif(X, P0, is_planar, n);
/*
 * Input:   A curve X,
 *          a point P0 on X,
 *          whether or not X is planar,
 *          and a precision n.
 * Output:  A development to precision n of P in a uniformizing parameter.
 *          The correct branch at P0 is chosen,
 *          and a coordinate is used in the case of a plane curve.
 */

if not is_planar then
    // TODO: Name issues.
    xs := GeneratorsSequence(CoordinateRing(Ambient(X)));
    K<x, y> := FunctionField(X);
    Pl0 := Place(P0);
    return [ Expand(K ! c, Pl0 : AbsPrec := n) : c in xs ];
end if;
f := DefiningEquations(X)[1];
R<x, y> := Parent(f);
xy := [x, y];
F := BaseRing(R);
RP<t> := PowerSeriesRing(F, n);
u := MyUniformizer(X, P0, is_planar);
S := PolynomialRing(RP);
/* The following is perhaps a bit ridiculous; copying and pasting would also have worked. */
if u eq x then
    i := 1;
else
    i := 2;
end if;
subst := [S!0, S!0];
subst[i] := t + P0[i];
subst[(i mod 2) + 1] := S.1;
h := hom<R -> S | subst>;
g := h(f);
for rt in [ tup[1] : tup in Roots(g) ] do
    if Coefficient(rt, 0) eq P0[(i mod 2) + i] then
        dev := [ RP!0, RP!0 ];
        dev[i] := t + P0[i];
        dev[(i mod 2) + 1] := rt;
        return dev;
    end if;
end for;

end function;


function NormalizeDiffBasis(X, P0, is_planar, B);
/*
 * Input:   A curve X,
 *          a point P0 on X,
 *          whether or not X is planar,
 *          and a basis of differentials B of X.
 * Output:  A differential basis Bnorm that is normalized with respect to the uniformizing parameter,
 *          the expansion of Bnorm,
 *          along with a matrix T such that multiplication by T sends B to Bnorm.
 */

g := #B;
P := DevelopInUnif(X, P0, is_planar, g);
du := Differential(MyUniformizer(X, P0, is_planar));
Q := FieldOfFractions(CoordinateRing(Ambient(X)));
BP := [ Evaluate(Q ! (b / du), P) : b in B ];
T := Matrix([ [ Coefficient(BP[i], j-1) : j in [1..g] ] : i in [1..g] ])^(-1);
Bnorm := [ &+[ T[i,j] * B[j] : j in [1..g] ] : i in [1..g] ];
return Bnorm, T;

end function;


function RamificationIndex(M);
/*
 * Input:   A matrix M that represents an endomorphism,
 *          after normalizing to an upper triangular form.
 * Output:  The ramification index of the corresponding Puiseux expansions.
 */

g := #Rows(M);
e := g;
while true do
    n := g div e;
    test := &or[ M[e*i, i] ne 0 : i in [1..n] ];
    if test then
        return e;
    end if;
    e := e - 1;
end while;

end function;


function InitializeImageBranch(M);
/*
 * Input:   A matrix M that represents an endomorphism,
 *          after normalizing to an upper triangular form,
 *          and a base point P0.
 * Output:  The leading coefficients of the corresponding Puiseux expansions.
 */

/* Recovering old invariants: */
F := Parent(M[1,1]);
g := #Rows(M);
e := RamificationIndex(M);

/* Normalized equations (depend only on the matrix): */
A := AffineSpace(F, g);
RA := CoordinateRing(A);
eqs := [ ];
for n in [1..g] do
    powersum := &+[ RA.i^n : i in [1..g] ];
    if n mod e eq 0 then
        Append(~eqs, powersum - e * M[n, n div e]);
    else
        Append(~eqs, powersum);
    end if;
end for;
S := Scheme(A, eqs);

/* Note: the upcoming steps are taken to avoid the use of an algebraic closure. */
RF<t> := PolynomialRing(F);
hc := [ RF!0 : i in [1..g] ];
hc[#hc] := t;
h := hom<RA -> RF | hc>;

/* By symmetry, this extension always suffices */
G := GroebnerBasis(ideal<RA | eqs>);
K := SplittingField(h(G[#G]));

/* Extending and evaluating: */
SK := BaseExtend(S, K);
Pt := Eltseq(Points(SK)[1]);
RK<tK> := PuiseuxSeriesRing(K, 1);
wK := tK^(1/e);
return [ Pt[i] * wK : i in [1..g] ];

end function;


function LinTermMult(lt1, lt2);
// Multiplies linear contributions

n := #Rows(Transpose(lt1));
return Matrix([ [ &+[ lt1[1, j] * lt2[1, i + 1 - j] : j in [1..i] ] : i in [1..n] ] ]);

end function;


function EvaluateOldGuess(fis, tj_old, n, m);
/*
 * Input:   A list of power series fis with relative precision at least 2n,
 *          a Puiseux series tj_old to be substituted in xj of relative precision n/e,
 *          and the integer n used above (for Magma problems with precision).
 * Output:  The evaluation of the first terms
 *          and the linear coefficients for the terms of higher order in the fi (tj).
 * NOTE:    The techniques are specific to the power series considered in the paper.
 */
/* TODO: Exploit symmetry when using this later? */

g := #fis;
e := Denominator(Valuation(tj_old));
F := BaseRing(Parent(tj_old));
RP<t> := PuiseuxSeriesRing(F, n + m);

/* Carrying over old guess to higher precision and extracting required information: */
tj := &+[ Coefficient(tj_old, i/e) * t^(i/e) : i in [1..n] ];
power0 := RP ! 1;
power1 := tj;
power_term0 := Matrix([ [ F ! 0 : k in [1..m] ] ]);
power_term1 := Matrix([ [ F ! 1 ] cat [ F ! 0 : k in [2..m] ] ]);
deriv := Matrix([ [ Coefficient(tj, k/e) : k in [1..m] ] ]);

/* Calculating contributions from higher powers of tj: */
powers := [ power0, power1 ];
power_terms := [ power_term0, power_term1 ];
power := power1;
power_term := power_term1;
for k in [2..((g - 1) + (n + m - 1))] do
    Append(~powers, power1 * powers[#powers]);
    if k le ((g - 1) + (m - 1)) then
        power_term := LinTermMult(deriv, power_term);
        Append(~power_terms, k * power_term);
    end if;
end for;

/* Combining for the given fis: */
evs0 := [ ];
lin_terms := [ ];
for i:=1 to g do
    ev0 := &+[ Coefficient(fis[i], k) * powers[k + 1] : k in [0..((i - 1) + (n + m - 1))] ];
    Append(~evs0, ev0);
    lin_term := Coefficient(fis[i], i - 1) * power_terms[i];
    /* Terms from higher powers: */
    for k:=2 to m do
        c := Coefficient(fis[i], (i - 1) + (k - 1));
        power_term := power_terms[i + k - 1];
        /* Add shift: */
        for l:=k to m do
            lin_term[1, l] +:= c * power_term[1, l + 1 - k];
        end for;
    end for;
    Append(~lin_terms, lin_term);
end for;

return evs0, lin_terms;

end function;


function MyEvaluate(X, P, is_planar, tj);
// Slightly faster evaluation in planar case.
// TODO: Similar to what went before, and we may wish to combine these.

if not is_planar then
    return [ Evaluate(c, tj) : c in P ];
end if;
f := DefiningEquations(X)[1];
R<x, y> := Parent(f);
xy := [x, y];
F := BaseRing(R);
RP<t> := Parent(tj);
S := PolynomialRing(RP);
P0 := [ Coefficient(c, 0) : c in P ];
u := MyUniformizer(X, P0, is_planar);
/* The following is perhaps a bit ridiculous; copying and pasting would also have worked. */
if u eq x then
    i := 1;
else
    i := 2;
end if;
subst := [S!0, S!0];
subst[i] := tj;
subst[(i mod 2) + 1] := S.1;
h := hom<R -> S | subst>;
g := h(f);
for rt in [ tup[1] : tup in Roots(g) ] do
    if Coefficient(rt, 0) eq Coefficient(P[(i mod 2) + i], 0) then
        dev := [ RP!0, RP!0 ];
        dev[i] := tj;
        dev[(i mod 2) + 1] := rt;
        return dev;
    end if;
end for;

end function;


function MySum(L);

if #L eq 0 then
    return 0;
end if;
return &+L;

end function;


function NthApproxs(X, P, is_planar, B, M, N);
/*
 * Input:   A curve X,
 *          a branch P,
 *          a basis of differentials B on X,
 *          a matrix representation M,
 *          and a precision N.
 * Output:  The Puiseux expansions of alpha (P) up to precision N.
 */

/* Recovering old invariants: */
g := #Rows(M);

/* Initializing and recovering corresponding Vandermonde matrix: */
ts0 := InitializeImageBranch(M);
e := Denominator(Valuation(ts0[1]));
A := Matrix([ [ Coefficient(ts0[j], 1/e)^(i - 1) : j in [1..g] ] : i in [1..g] ])^(-1);
F := BaseRing(Parent(ts0[1]));

/* Evaluating differentials to power series: */
P0 := X ! [ Coefficient(c, 0) : c in P ];
du := Differential(MyUniformizer(X, P0, is_planar));
Q := FieldOfFractions(CoordinateRing(Ambient(X)));
BP := [ Evaluate(Q ! (b / du), P) : b in B ];
ev_rhs := [ &+[ M[i, j] * BP[j] : j in [1..g] ] : i in [1..g] ];

/* Lifting: */
ts := ts0;
n := 1;
while n lt N do
    m := Minimum(n, N - n);
    // TODO: Clinch relative precision here.
    R<t> := PuiseuxSeriesRing(F, n + m);
    evss0 := [ ];
    lin_termss := [ ];
    for tj in ts do
        evs0, lin_terms := EvaluateOldGuess(BP, tj, n, m);
        Append(~evss0, [ R ! ev0 : ev0 in evs0 ]);
        Append(~lin_termss, lin_terms);
    end for;
    ts := [ &+[ Coefficient(ts[i], j/e) * t^(j/e) : j in [1..n] ] : i in [1..g] ];
    ev_lhs := [ &+[ evss0[j][i] * Derivative(ts[j]) : j in [1..g] ] : i in [1..g] ];
    diffs := [ ev_rhs[i] - ev_lhs[i] : i in [1..g] ];
    /* Yes, this is a horrible stack of for-loops. But it is fast. */
    for k in [1..m] do
        diffsdiffed := [ Coefficient(diffs[i], (i - 1 + n + k)/e - 1) : i in [1..g] ];
        for i in [1..g] do
            for j in [1..g] do
                for l in [1..k] do
                    /* Next line has zero entry at k: */
                    diffsdiffed[i] -:= ((n + l)/e) * Coefficient(ts[j], (n + l)/e) * Coefficient(evss0[j][i], (i - 1 + k - l)/e);
                    gamma := MySum([ lin_termss[j][i][1, l + 1 - ind] * Coefficient(ts[j], (n + ind)/e) : ind in [1..Minimum(l, k - 1)] ]);
                    diffsdiffed[i] -:= ((k + 1 - l)/e) * Coefficient(ts[j], (k + 1 - l)/e) * gamma;
                end for;
            end for;
        end for;
        v := Matrix([ [ (e/(i - 1 + n + k)) * diffsdiffed[i] ] : i in [1..g] ]);
        w := A*v;
        ts := [ ts[i] + w[i, 1] * t^((n + k)/e) : i in [1..g] ];
    end for;
    n *:= 2;
end while;

R<t> := PuiseuxSeriesRing(F, n + m + g - 1);
return [ MyEvaluate(X, P, is_planar, R ! t) : t in ts ];

end function;


function NthApproxsOld(X, P, is_planar, B, M, N);
/*
 * Input:   A curve X,
 *          a branch P,
 *          a basis of differentials B on X,
 *          a matrix representation M,
 *          and a precision N.
 * Output:  The Puiseux expansions of alpha (P) up to precision N.
 */

/* Recovering old invariants: */
g := #Rows(M);

/* Initializing and recovering corresponding Vandermonde matrix: */
ts0 := InitializeImageBranch(M);
e := Denominator(Valuation(ts0[1]));
A := Matrix([ [ Coefficient(ts0[j], 1/e)^(i - 1) : j in [1..g] ] : i in [1..g] ])^(-1);
F := BaseRing(Parent(ts0[1]));

/* Evaluating differentials to power series: */
P0 := X ! [ Coefficient(c, 0) : c in P ];
du := Differential(MyUniformizer(X, P0, is_planar));
Q := FieldOfFractions(CoordinateRing(Ambient(X)));
BP := [ Evaluate(Q ! (b / du), P) : b in B ];

/* Lifting: */
ts := ts0;
n := 1;
while n lt N do
    n := n + 1;
    // TODO: Check that this is enough precision in genus > 2.
    // (Issue is whether precision gets pushed out properly.)
    RP<tP> := PuiseuxSeriesRing(F, n);
    wP := tP^(1/e);
    ts := [ &+[ Coefficient(ts[i], j/e) * wP^j : j in [1..n-1] ] : i in [1..g] ];
    ev_lhs := [ &+[ Evaluate(BP[i], ts[j]) * Derivative(ts[j]) : j in [1..g] ] : i in [1..g] ];
    ev_rhs := [ &+[ M[i,j] * BP[j] : j in [1..g] ] : i in [1..g] ];
    diffs := [ ev_rhs[i] - ev_lhs[i] : i in [1..g] ];
    v := Matrix([ [ (e / (n + i - 1)) * Coefficient(diffs[i], (n + i - 1)/e - 1) ] : i in [1..g] ]);
    w := A*v;
    ts := [ ts[i] + w[i, 1] * tP^(n/e) : i in [1..g] ];
end while;

R<t> := PuiseuxSeriesRing(F, n + g - 1);
return [ MyEvaluate(X, P, is_planar, R ! t) : t in ts ];

end function;


function MyTest();

n := 4;
m := 4;
g := 3;
// TODO: Clinch relative precision here.
RP<x> := PowerSeriesRing(Rationals(), 2*n + g - 1);
fi1 := &+[ l*x^(l-1) : l in [1..2*n] ];
fi2 := &+[ l^2*x^(l-1) : l in [2..2*n] ];
fi3 := &+[ l^3*x^(l-1) : l in [3..2*n] ];
fis := [ fi1, fi2, fi3 ];
K<a,b,c,d> := RationalFunctionField(Rationals(), n);
RP<t> := PuiseuxSeriesRing(K, 2*n);
tj := 111*t^(1/2) + 17*t^(2/2) - 11*t^(3/2) + 31*t^(4/2) + a*t^(5/2) + b*t^(6/2) + c*t^(7/2) + d*t^(8/2);
ev0, M := EvaluateOldGuess(fis, tj, n, m);

return [* [ tj^k : k in [1..2*n] ], [ Evaluate(fi, tj) : fi in fis ], ev0, M *];

end function;
