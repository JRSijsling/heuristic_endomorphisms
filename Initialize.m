/***
 *  Initialization of the Magma part of the package
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */

// Global constants; should perhaps just copy these
prec0 := 200;
epscomp0 := 10^(-prec0 + 30);
epsLLL0 := 5^(-prec0 + 7);
epsinv0 := 2^(-4);

// Linear algebra routines
forward NumericalSolve;
forward InvertibleSubmatrix;
forward SplitPeriodMatrix;
forward CombinePeriodMatrix;
forward IntegralKernel;
//forward OptimizedPeriodMatrix;
forward ConjugateMatrix;
forward MatrixInBasis;
forward MatrixRatInBasisOverNF;

// Functionality for recognizing complex numbers as algebraic numbers and
// related optimizations
forward CompositeReduce;
forward MySplittingField;
forward PolynomializeElement;
forward PolynomializeMatrix;
forward FractionalApproximation;
forward AlgebraizeElementInField;
forward NearbyRoot;
forward AlgebraizeMatricesInField;
forward IntegralRepresentationNF;
//forward PartialLMFDBLabel;

// Finding an approximate basis for the geometric endomorphism ring through LLL
forward ComplexStructure;
forward RationalEndomorphismEquations;
forward AnalyticRepresentation;
forward PeriodMatrix;
forward GeometricEndomorphismBasisFromPeriodMatrix;

// Algebraizing the basis
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
forward SatoTateGroupG2;
forward DecompositionIdempotentsG2;
forward GeoEndRRShorthand;

// Determining the Sato-Tate group
forward SatoTateGroupG2;

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

// The actual loading
AttachSpec("spec");
load "heuristic/Linear.m";
load "heuristic/Recognition.m";
load "heuristic/Analytic.m";
load "heuristic/Algebraic.m";
load "heuristic/SatoTate.m";
load "heuristic/Canonize.m";
load "heuristic/Decomposition.m";
load "polarization/Rosati.m";

load "heuristic/AlgebraicOneOff.m";

