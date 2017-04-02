/***
 *  Linear algebra routines
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


function NumericalLeftSolve(A, B);
// Input:   An invertible matrix A and a matrix B of compatible dimensions.
// Output:  An approximation of the matrix X such that X * A = B.

// An alternative that should be discussed is to find a full basis at once by
// stacking all of -B under A. This leads to inversion though.

// Yet another option is to do this for matrices of full rank A by piling zero
// matrices on top of the stack above and using the echelonized return form of
// the kernel. This does not work well with NumericalKernel though.

// FIXME: It seems that this causes instabilities sometimes; I therefore revert
// to naive inversion.

return B * A^(-1);

m := #Rows(Transpose(B));
XRows := [];
for r in Rows(B) do
    Br := Matrix([r]);
    Br := VerticalJoin(A, -Br);
    // Adding a trivial column to make NumericalKernel work;
    // the location of this column at the beginning is actually important!
    Br := HorizontalJoin(ZeroMatrix(BaseRing(A), m + 1, 1), Br);
    K := NumericalKernel(Br);
    v := K[1];
    Xr := Matrix([[ v[i] / v[m + 1] : i in [1..m] ]]);
    Append(~XRows, Xr);
end for;

X := VerticalJoin(XRows);
return X;

end function;


function InvertibleSubmatrix(P : epsinv := epsinv0);
// Input:   A period matrix P of dimension 2g x g along with a tolerance epsinv
//          for testing triviality of subdeterminants.
// Output:  An invertible submatrix and the corresponding subset of {1..2g}. We
//          do not attempt to find the best one.

g := #Rows(Transpose(P));
RP := Rows(P);
for s in Subsets({1..2*g}, g) do
    s0 := s;
    P0 := Matrix([RP[i] : i in s]);
    if Abs(Determinant(P0)) gt epsinv then
        return P0, s0;
    end if;
end for;

error Error("Failed to find an invertible submatrix");

end function;



function SubmatrixOfRank(P, TargetRank, ColumnsOrRows : epsinv := epsinv0);
// Input:   A period matrix P of dimension r x c along with a tolerance epsinv
//          for testing triviality of subdeterminants.
//          The number of rows/columns to extract, whether to extract rows or
//		 columns.
// Output:  A TargetRank x #columns submatrix M and the corresponding subset 
//          of {1..#rows/#columns}. The matrix M should have rank TargetRank
//          We do not attempt to find the best one.


Flipped := false;
if ColumnsOrRows eq "Columns" then
    Flipped := true;
    P := Transpose(P);
end if;

r := #Rows(P);
c := #Rows(Transpose(P));


RP := Rows(P);
for s in Subsets({1..r}, TargetRank) do
    s0 := s;
    P0 := Matrix([RP[i] : i in s]);
    if NumericalRank(P0 : Epsilon:=RealField(100)!epsinv) eq TargetRank then 
        if not Flipped then
            return P0, s0;
        else
            return Transpose(P0), s0;
        end if;
    end if;
end for;

error Error("Failed to find submatrix with the desired rank");

end function;


function SplitPeriodMatrix(P);
// Input:   A period matrix P of dimension 2g x g.
// Output:  Its real and imaginary parts stacked on top of each other in a 2g x
//          2g matrix.

PSplitRe := Matrix([[Real(c) : c in Eltseq(r)] : r in Rows(P)]);
PSplitIm := Matrix([[Im(c) : c in Eltseq(r)] : r in Rows(P)]);
PSplit := HorizontalJoin([PSplitRe, PSplitIm]);

return PSplit;

end function;


function CombinePeriodMatrix(PSplit);
// Input:   A split period matrix as in the output of the previous algorithm.
// Output:  The corresponding original period matrix.

// Basic invariants
c := #Rows(Transpose(PSplit));
r := #Rows(PSplit);

k0 := BaseRing(PSplit);
k<i> := ComplexField(k0);

// Combining parts
PRe := Matrix(k, Submatrix(PSplit, [1..r], [1..(c div 2)]));
PIm := Matrix(k, Submatrix(PSplit, [1..r], [(c div 2)+1..c]));
P := PRe + i*PIm;

return P;

end function;


function IntegralLeftKernel(M : epsLLL := epsLLL0);
// Input:   A real matrix whose rows are to be cancelled by integral
//          combinations and an LLL parameter epsLLL.
// Output:  The corresponding LLL approximation of the integral left kernel.

// Usual LLL approximation trick
MI := IdentityMatrix(BaseRing(M), #Rows(M));
Ma := HorizontalJoin(MI, (1 / epsLLL) * M);
L, K := LLL(Ma);

return K;

end function;


function ConjugateMatrix(sigma, M);
// Input:   A matrix M over a field and an automorphism sigma of that field.
// Output:  The matrix M after conjugating its coefficients by sigma.

return Matrix([[sigma(elt) : elt in Eltseq(row)] : row in Rows(M)]);

end function;


function MatrixInBasis(M, Bs);
// Input:   A matrix M in a matrix algebra of which Bs form a basis
// Output:  A representation of M in that basis

MM := Matrix([Eltseq(M)]);
MBs := Matrix([Eltseq(B) : B in Bs]);
V := Solution(MBs, MM);

return Matrix(V);

end function;


function MatrixRatInBasisOverNF(M, Bs);
// Input:   A matrix M in a matrix algebra of which Bs form a basis
// Output:  A representation of M in that basis

MM := Matrix([&cat[Eltseq(m) : m in Eltseq(M)]]);
MBs := Matrix([&cat[Eltseq(b) : b in Eltseq(B)] : B in Bs]);
V := Solution(MBs, MM);

return Matrix(V);

end function;


function SaturateLattice(FullLattice, SubLattice : epscomp := epscomp, epsLLL := epsLLL)
// Input:   Matrices FullLattice, SubLattice whose rows generate lattices
//		 such that the latter has finite index in the former
// Output:  A matrix of the same size as SubLattice whose rows generate the
//		 same lattice as FullLattice
    MatrixSize := Min(#Rows(SubLattice), #Rows(Transpose(FullLattice)));
    
    SaturationMatrix, SatColumns := SubmatrixOfRank(SubLattice, MatrixSize, "Columns");
    LMat := Matrix(Rationals(), 0,MatrixSize, []);
    for row in Rows(FullLattice) do
        linearDependence := NumericalLeftSolve(SaturationMatrix, Matrix([[row[j] : j in SatColumns]]) );
        linearDependence := Matrix([ [ FractionalApproximation(c : epscomp := epscomp, epsLLL := epsLLL) : c in Eltseq(row) ] : row in Rows(linearDependence) ]);
        LMat := VerticalJoin(LMat, linearDependence);
    end for;
    
    M := Matrix(Basis(Lattice(LMat)));
    return Matrix(BaseRing(FullLattice),M)*SubLattice;
end function;
