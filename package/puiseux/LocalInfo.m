forward LiftPuiseuxSeries;

forward InitializeMatrix;
forward PuiseuxRamificationIndex;
forward InitializeImageBranch;

forward RootWithHensel;
forward DevelopPoint;

forward InitializeLift;
forward CreateLiftIterator;
forward ApproximationsFromTangentAction;


function LiftPuiseuxSeries(f, PR, e)

L := [ Coefficient(f, i/e)*PR.1^(i/e) : i in [0..e*AbsolutePrecision(f) - 1] ];
if #L eq 0 then
    return PR ! 0;
end if;
return PR ! (&+L);

end function;


function PuiseuxRamificationIndex(M)
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


function InitializeImageBranch(M)
/*
 * Input:   A matrix M that represents an endomorphism
 *          after normalizing to an upper triangular form.
 * Output:  The leading coefficients of the corresponding Puiseux expansions
 *          and the polynomial that gives their field of definition.
 */

/* Recovering old invariants: */
F := Parent(M[1,1]);
g := #Rows(M);
e := PuiseuxRamificationIndex(M);

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

/* The upcoming steps are taken to avoid the use of an algebraic closure */
RF := PolynomialRing(F);
hc := [ RF!0 : i in [1..g] ];
hc[#hc] := RF.1;
h := hom<RA -> RF | hc>;

/* By symmetry, this extension always suffices */
G := GroebnerBasis(ideal<RA | eqs>);
if not IsFinite(F) then
    K := GaloisSplittingField(h(G[#G]));
else
    K := SplittingField(h(G[#G]));
end if;

/* Extending and evaluating: */
SK := BaseExtend(S, K);
P := Eltseq(Points(SK)[1]);
if e eq 1 then
    RK := PowerSeriesRing(K, 2);
    wK := RK.1;
else
    RK := PuiseuxSeriesRing(K, 2);
    wK := RK.1^(1/e);
end if;
return [ P[i] * wK + O(wK^2) : i in [1..g] ], h(G[#G]);

end function;


function RootWithHensel(f, P0, n)
/*
 * Input:   An algebraic relation f between two variables,
 *          a point P0 that satisfies this,
 *          and the requested number of digits n.
 * Output:  A corresponding development of both components.
 *
 * The relation f has to be non-singular when developing in y,
 * and the x-coordinate can be specified as a Puiseux series. */

if Type(P0[1]) in [ RngSerPuisElt, RngSerPowElt ] then
    F := BaseRing(Parent(P0[1]));
else
    F := Parent(P0[1]);
end if;
df := Derivative(f, 2);
x0 := P0[1]; y0 := P0[2];

PR := PuiseuxSeriesRing(F, 1);
x := PR ! Coefficient(PR ! x0, 0); y := PR ! Coefficient(PR ! y0, 0);
if n eq 0 then
    return [x, y];
end if;

log := 0;
while log le Ceiling(Log(2, n)) - 1 do
    prec := Minimum(2^(log + 1), n);
    PR := PuiseuxSeriesRing(F, prec);
    if x0 in F then
        e := 1;
        x := x0 + PR.1;
    else
        /* TODO: We further assume a somewhat special kind of series here,
         *       which we encounter in our applications. */
        e := Integers() ! (1/Valuation(x0 - Coefficient(x0, 0)));
        x := LiftPuiseuxSeries(x0, PR, e);
    end if;
    y := LiftPuiseuxSeries(y, PR, e);
    h := -Evaluate(f, [x, y])/Evaluate(df, [x, y]);
    y +:= h;
    log +:= 1;
end while;
return [x, y];

end function;


function DevelopPoint(X, P0, n)
/*
 * Input:   A curve X,
 *          a point P0 on it (may be over a series ring),
 *          and a precision n.
 * Output:  A development to precision n of P in a uniformizing parameter.
 *          The correct branch at P0 is chosen,
 *          and a coordinate is used in the case of a plane curve.
 */

if not X`is_planar then
    /* Here only for constant points, in which case we fall back to the given
     * base point. We do get an expansion that may not be in our uniformizer. */
    return [ Expand(X`K ! c, Place(X`Q0) : AbsPrec := n) : c in GeneratorsSequence(X`R) ];
end if;
f := X`DEs[1];
if X`index eq 1 then
    return RootWithHensel(f, P0, n);
else
    f_swap := Evaluate(f, [(X`R).2, (X`R).1]);
    P0_swap := [ P0[2], P0[1] ];
    P := RootWithHensel(f_swap, P0_swap, n);
    return [ P[2], P[1] ];
end if;

end function;


function InitializeLift(X, M)

P0 := X`Q0;
e := PuiseuxRamificationIndex(M);
tjs0 := InitializeImageBranch(M);
PR := Parent(tjs0[1]);

/* Creating P */
P := [ PR ! P0[1], PR ! P0[2] ];
P[X`index] +:= PR.1;
Pnew := [ PR ! c : c in DevelopPoint(X, P, 2) ];
P := [ PR ! c : c in DevelopPoint(X, P, 2) ];

/* Creating Qjs */
Qjs := [ [ PR ! P0[1], PR ! P0[2] ] : i in [1..X`g] ];
for i in [1..X`g] do
    Qjs[i][X`index] +:= tjs0[i];
end for;
Qjs := [ [ PR ! c : c in DevelopPoint(X, Qj, 2) ] : Qj in Qjs ];
return P, Qjs;

end function;


function CreateLiftIterator(X, M)

index_unif := X`index;
index_other := (index_unif mod 2) + 1;
f := X`DEs[1];
df := Derivative(f, index_other);
B := X`NormB;
g := X`g;

e := PuiseuxRamificationIndex(M);

    function Iterate(P, Qjs, n);

    /* Create ring of higher precision: */
    K := BaseRing(Parent(P[1]));
    prec := Minimum(2*Precision(Parent(P[1])), n);
    PR := PuiseuxSeriesRing(K, prec);

    /* Lift Qjs: */
    Qjs := [ [ LiftPuiseuxSeries(c, PR, e) : c in Qj ] : Qj in Qjs ];
    for i in [1..g] do
        h := -Evaluate(f, Qjs[i])/Evaluate(df, Qjs[i]);
        Qjs[i][index_other] +:= h;
    end for;

    /* Calculate P to higher precision: */
    P := [ LiftPuiseuxSeries(c, PR, 1) : c in P ];
    P[index_unif] := Coefficient(P[index_unif], 0) + PR.1;
    h := -Evaluate(f, P)/Evaluate(df, P);
    P[index_other] +:= h;

    /* Calculate LHS: */
    dtjs := [ Derivative(Qj[index_unif]) : Qj in Qjs ];
    BPijs := Matrix([ [ Evaluate(B[i], Qjs[j]) : j in [1..g] ] : i in [1..g] ]);
    F_ev := Matrix([ [ Integral(&+[ BPijs[i,j] * dtjs[j] : j in [1..g] ]) : i in [1..g] ] ]);
    DF_ev := Transpose(BPijs);

    /* Calculate RHS: */
    BQjs := [ Evaluate(B[j], P) : j in [1..g] ];
    G_ev := Matrix(PR, [ [ Integral(&+[ M[i,j] * BQjs[j] : j in [1..g] ]) : i in [1..g] ] ]);

    /* Calculate Hensel correction: */
    H := -(F_ev - G_ev) * DF_ev^(-1);

    /* Calculate Qjs to higher precision: */
    for i in [1..g] do
        Qjs[i][index_unif] +:= H[1,i];
        h := -Evaluate(f, [Qjs[i][1], Qjs[i][2]])/Evaluate(df, [Qjs[i][1], Qjs[i][2]]);
        Qjs[i][index_other] +:= h;
    end for;
    return P, Qjs;

    end function;

return Iterate;

end function;


intrinsic ApproximationsFromTangentAction(X::Crv, M::AlgMatElt, n::RngIntElt) -> Tup, Tup
{Given a curve X, a matrix M that gives the tangent representation of an
endomorphism on the normalized basis of differentials, and an integer n,
returns a development of the branches to precision at least O(n).}

e := PuiseuxRamificationIndex(M);
P, Qjs := InitializeLift(X, M);
IterateLift := CreateLiftIterator(X, M);
/* Next line avoids precision loss */
/* TODO: Determine exact bound needed here */
for i:= 1 to Ceiling(Log(2, n + e + 1)) do
    P, Qjs := IterateLift(P, Qjs, n + e + 1);
end for;
return P, Qjs;

end intrinsic;

