function SplittingInfoOneOff(Rs)

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
return Ds;

end function;

