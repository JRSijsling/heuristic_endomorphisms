/***
 *  Elliptic curves from decompositions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic DecompositionDenominator(R::.) -> RngIntElt
{Degree of a projection morphism.}

return LCM([ Denominator(c) : c in Eltseq(R) ]);

end intrinsic;


intrinsic Idempotents(EndList::List) -> SeqEnum
{Idempotents of the endomorphism algebra.}

return 0;

end intrinsic;


function SmallIdempotent(Q, GensQ);
// Input:   A split quaternion algebra and a basis for an order.
// Output:  A nilpotent element that is small with respect to that order in the
//          sense that the numerator is small

test, M, f := IsMatrixRing(Q : Isomorphism := true);
invf := Inverse(f);

// TODO: What follows is obviously stupid (it does not even use GensQ) and
// should be improved.
return invf(M![1,0,0,0]);

end function;


intrinsic AnalyticProjection(idem::., P::.) -> List
{From an idempotent, extracts corresponding lattice and projection.}
    // Extract the complex field
    CC := BaseRing(P); RR := RealField(CC);

    // Create analytic idempotent and project:
    PEllHuge := P * idem;

    // Compute rank of projection
    gQuotient := NumericalRank(PEllHuge : Epsilon := RR`epsinv);

    PEllHugeSplit := SplitPeriodMatrix(PEllHuge);

    PEllBigSplit, col_numbers := SubmatrixOfRank(PEllHugeSplit, 2*gQuotient, "Rows"); // extract 2g independent rows
    PEllBigSplit := SaturateLattice(PEllHugeSplit, PEllBigSplit : epscomp := epscomp, epsLLL := epsLLL); // ensure that these rows generate the full lattice

    PEllBig := CombinePeriodMatrix(PEllBigSplit); // go back to the complex representation

    PreliminaryLatticeMatrix := SubmatrixOfRank(PEllBig, gQuotient, "Columns" : epsinv := epsinv); // extract g columns (i.e. decide which projection to use)

    PreliminaryLatticeMatrix := Transpose(PreliminaryLatticeMatrix); // necessary before calling SaturateLattice
    PEllBig := Transpose(PEllBig);
    LatticeMatrix := Transpose(SaturateLattice(PEllBig, PreliminaryLatticeMatrix : epscomp := epscomp, epsLLL := epsLLL));
end intrinsic;


intrinsic AnalyticProjections(idems::., P::.) -> List
{From idempotents, extracts corresponding lattices and projections.}

return [* AnalyticProjection(idem, P) : idem in idems *];

end intrinsic;
