function SplittingInfoOneOff(AsAlg, As, Rs)

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

Ds, idemsC := DirectSumDecomposition(C);
idemsC := [ [ Rationals() ! c : c in Eltseq(idem) ] : idem in idemsC ];
GensC := [ [ Rationals() ! c : c in Eltseq(gen) ] : gen in GensC ];
idemsRep := [ Eltseq(MatrixInBasis(idem, GensC)) : idem in idemsC ];
N := #AsAlg;
idemsAlg := [ &+[ idem[n] * AsAlg[n] : n in [1..N] ] : idem in idemsRep ];
idemsAn := [ &+[ idem[n] * As[n] : n in [1..N] ] : idem in idemsRep ];
idemsRat := [ &+[ idem[n] * Rs[n] : n in [1..N] ] : idem in idemsRep ];

return Ds, idemsAlg, idemsAn;

end function;

