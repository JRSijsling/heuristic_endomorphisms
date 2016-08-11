/***
 *  Initialization of the Magma part of the package
 *
 *  Copyright (C) 2016  J.R. Sijsling (sijsling@gmail.com)
 *
 *  Distributed under the terms of the GNU General License (GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc., 51
 *  Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
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
forward CanonizeMatrix;

// Elliptic curves from decompositions
forward DecompositionDegree;
forward SmallIdempotent;
forward InducedEmbedding;
forward LatticesFromIdempotents;

// Explicit projection morphisms
forward ProjectionToEllipticFactorG2;

// Rosati involution functionality
forward RosatiInvolution;

// The actual loading
load "Linear.m";
load "Recognition.m";
load "Analytic.m";
load "Algebraic.m";
load "SatoTate.m";
load "Canonize.m";
load "Decomposition.m";
load "Projection.m";
load "Rosati.m";
