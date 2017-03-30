/***
 *  Initialization of the Magma part of the package
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

// NOTE: The end user should only overwrite these options in Wrapper.sage
prec0 := 300;
epscomp0 := 10^(-prec0 + 30);
epsLLL0 := 5^(-prec0 + 7);
epsinv0 := 2^(-prec0 + 30);

// Relative number fields
forward EmbedAsComplexPolynomials;
forward EmbedAsComplexPolynomial;

// Periods
forward PeriodMatrixHyperelliptic;
forward PeriodMatrixPlane;

// Linear algebra routines
forward NumericalLeftSolve;
forward InvertibleSubmatrix;
forward SplitPeriodMatrix;
forward CombinePeriodMatrix;
forward IntegralKernel;
forward ConjugateMatrix;
forward MatrixInBasis;
forward MatrixRatInBasisOverNF;

// Functionality for recognizing complex numbers as algebraic numbers and
// related optimizations
// TODO: Quite ugly right now
forward CompositeReduce;
forward MySplittingField;
forward PolynomializeElement;
forward PolynomializeMatrix;
forward FractionalApproximation;
forward AlgebraizeElementInField;
forward NearbyRoot;
forward AlgebraizeMatricesInField;
forward IntegralRepresentationNF;

// Finding an approximate basis for the geometric endomorphism ring through LLL
forward ComplexStructure;
forward RationalEndomorphismEquations;
forward AnalyticRepresentation;
forward GeometricEndomorphismBasisFromPeriodMatrix;

// TODO: Generalize and modularize
forward EndomorphismBasisOverSubfield;
forward CompareSubfield;
forward EndomorphismBasisOverField;

forward EndomorphismLatticeG2SingleElement;
forward EndomorphismLatticeG2;

forward EndomorphismData;
forward EndomorphismAlgebraFactorsGeneric;
forward EndomorphismRingGeneric;
forward EndomorphismAlgebraFactorsQQG2;
forward EndomorphismAlgebraFactorsRRG2;
forward EndomorphismRingG2;
forward DecompositionIdempotentsG2;
forward GeoEndRRShorthand;

// Determining the Sato-Tate group in genus 2
forward SatoTateGroupG2;

// To genus 3
forward EndomorphismLatticeG3SingleElement;
forward EndomorphismLatticeG3;

forward EndomorphismDataG3;
forward EndomorphismAlgebraFactorsQQG3;
forward EndomorphismAlgebraFactorsRRG3;
forward EndomorphismRingG3;

// Subfield canonization
forward CanonizeMatrices;

// Elliptic curves from decompositions
forward DecompositionDegree;
forward SmallIdempotent;
forward InducedEmbedding;
forward LatticesFromIdempotents;

// Rosati involution functionality
forward RosatiInvolution;

// Explicit projection morphisms
forward ProjectionToEllipticFactorG2;

// Heuristic part
AttachSpec("spec");
load "heuristic/Relative.m";
load "heuristic/Periods.m";
load "heuristic/Linear.m";
load "heuristic/Analytic.m";
load "heuristic/Recognition.m";
load "heuristic/OverField.m";
load "heuristic/Lattice.m";
load "heuristic/Data.m";
load "heuristic/SatoTate.m";
load "heuristic/LatticeG3.m";
load "heuristic/DataG3.m";
load "heuristic/Canonize.m";
load "heuristic/Decomposition.m";
load "heuristic/FromBig.m";

// Verification part
load "polarization/Rosati.m";
