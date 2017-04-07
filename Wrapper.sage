"""
 *  Class wrappers for the algorithms
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

def TypeTest(X):
    str0 = str(X.__class__)
    if str0 == "<class 'sage.schemes.curves.projective_curve.ProjectivePlaneCurve_with_category'>":
        return "plane"
    elif str0 == "<class 'sage.schemes.hyperelliptic_curves.hyperelliptic_g2_rational_field.HyperellipticCurve_g2_rational_field_with_category'>":
        return "hyperelliptic"
    elif str0 == "<class 'sage.schemes.hyperelliptic_curves.hyperelliptic_rational_field.HyperellipticCurve_rational_field_with_category'>":
        return "hyperelliptic"
    elif str0 == "<class 'sage.schemes.hyperelliptic_curves.hyperelliptic_g2_generic.HyperellipticCurve_g2_generic_with_category'>":
        return "hyperelliptic"
    else:
        return "generic"

class EndomorphismData:
    def __init__(self, X, prec, bound = 0, have_oldenburg = False):
        self.curve_type = TypeTest(X)
        if self.curve_type == "hyperelliptic":
            f, h = X.hyperelliptic_polynomials()
            self.f = magma(f)
            self.h = magma(h)
            self.base_field = magma.BaseRing(self.f)
            embedded_list = magma.EmbedAsComplexPolynomials([self.f, self.h], prec)
            self._fCC_, self._hCC_ = embedded_list
        elif self.curve_type == "plane":
            self.F = magma(X.defining_polynomial())
            self.base_field = magma.BaseRing(self.F)
            embedded_list = magma.EmbedAsComplexPolynomials([self.F], prec)
            self._FCC_ = embedded_list[1]
        self.prec = prec
        self.bound = bound
        self.have_oldenburg = have_oldenburg

    def __repr__(self):
        return ReprEndomorphismData(self)

    def period_matrix(self):
        if not hasattr(self, "_P_"):
            self._P_ = magma.PeriodMatrix(self._fCC_, self._hCC_, HaveOldenburg = self.have_oldenburg)
        return self._P_

    def geometric_representations(self):
        if not hasattr(self, "_geo_reps_"):
            self._P_ = self.period_matrix()
            self._geo_reps_ = magma.GeometricEndomorphismBasisApproximations(self._P_)
            self._AsPol_ = magma.RelativeMinimalPolynomialsMatrices(self._geo_reps_[1], self.base_field)
            self._endo_fod_ = Relative_Splitting_Field(self._AsPol_, bound = self.bound)
            self._geo_reps_ = magma.GeometricEndomorphismBasis(self._geo_reps_, self._endo_fod_)
            self._AsAlg_ = self._geo_reps_[1]
            self._Rs_ = self._geo_reps_[2]
            self._As_ = self._geo_reps_[3]
        return self._geo_reps_

    def endomorphism_field(self):
        self._geo_reps_ = self.geometric_representations()
        return self._endo_fod_

    def geometric(self):
        self._geo_reps_ = self.geometric_representations()
        return OverField(self, K = "geometric")

    def over_base(self):
        self._geo_reps_ = self.geometric_representations()
        return OverField(self, K = "base")

    def over_field(self, K):
        self._geo_reps_ = self.geometric_representations()
        return OverField(self, K = K)

    def lattice(self):
        if not hasattr(self, "_lat_"):
            self._geo_reps_ = self.geometric_representations()
            self._lat_ = magma.EndomorphismLattice(self._geo_reps_)
        return Lattice(self._lat_)

    def rosati_involution(self, A):
        self._geo_reps_ = self.geometric_representations()
        return magma.RosatiInvolution(self._geo_reps_, A)

    def degree_estimate(self, A):
        self._geo_reps_ = self.geometric_representations()
        return magma.DegreeEstimate(self._geo_reps_, A)

    def verify_algebra(self):
        return True

    def verify_saturated(self):
        self._geo_reps_ = self.geometric_representations()
        return magma.VerifySaturated(self._geo_reps_, self._P_)

    def verify_representations(self):
        self._geo_reps_ = self.geometric_representations()
        XL, P0L, AsL = magma.NonWeierstrassBasePointHyp(self.X, self.base_field, self._AsAlg_, nvals = 3)
        d = self.degree_estimate(AL)
        return magma.VerifyRepresentations(XL, P0L, ALs, LowerBound = 2*d + 2)

    def verify(self):
        return (self.verify_algebra() and self.verify_saturated() and self.verify_representations())

    def decomposition(self):
        self._P_ = self.period_matrix()
        self._lat_ = self.lattice()
        return Decomposition(self)

class Lattice:
    def __init__(self, lat):
        self._lat_ = lat
        self._lat_reps_ = lat[1]
        self._lat_structs_ = lat[2]
        self._lat_descs_ = lat[3]

    def __repr__(self):
        return ReprLattice(self)

    def representations(self):
        # TODO: Make what follows a conversion function, Galois group
        return self._lat_reps_

    def structures(self):
        # TODO: Make what follows a conversion function, Galois group
        return self._lat_structs_

    def descriptions(self):
        # TODO: Make what follows a conversion function, Galois group
        return self._lat_descs_

class OverField:
    def __init__(self, EndJac, K = "geometric"):
        self.X = EndJac.X
        self.base_field = EndJac.base_field
        self.field = K
        self._geo_reps_ = EndJac._geo_reps_

    def __repr__(self):
        return ReprOverField(self)

    def representations(self):
        if not hasattr(self, "_reps_"):
            if self.field == "geometric":
                self._reps_ = self._geo_reps_
            elif self.field == "base":
                self._reps_ = magma.EndomorphismBasis(self._geo_reps_, self.base_field)
            else:
                self._reps_ = magma.EndomorphismBasis(self._geo_reps_, magma(self.field))
        return self._reps_

    def representations_tangent(self):
        self._reps_ = self.representations()
        return self._reps_[1]

    def representations_homology(self):
        self._reps_ = self.representations()
        return self._reps_[2]

    def representations_approximations(self):
        self._reps_ = self.representations()
        return self._reps_[3]

    def structure(self):
        if not hasattr(self, "_struct_"):
            self._reps_ = self.representations()
            self._struct_, self._desc_ = magma.EndomorphismStructure(self._reps_, nvals = 2)
        return self._struct_

    def description(self):
        if not hasattr(self, "_desc_"):
            self._reps_ = self.representations()
            self._struct_, self._desc_ = magma.EndomorphismStructure(self._reps_, nvals = 2)
        return ReprDescription(self._desc_)

class Decomposition:
    # TODO: We could have different fields here. Do we want a catch-all?
    def __init__(self, EndJac):
        self.X = EndJac.X

    def __repr__(self):
        return ReprDecomposition(self)

    def field_of_definition(self):
        if not hasattr(self, "_decomp_fod_"):
            return 0
        # TODO: Run through lattice
        return 0

    def idempotents(self):
        if not hasattr(self, "_idems_"):
            return 0
        # TODO: Run through lattice; count over closure and go down, possibly even finer
        return 0

    def putative_factors(self):
        # TODO: Defer to Magma. Combination of two functions
        if not hasattr(self, "_factors_"):
            self._decomp_fod_ = self.field_of_definition()
            self._idems_ = self.idempotents()
            self._factors_ = magma.PutativeFactorsFromIdempotents(self._idems_)
        return self._factors_

    def morphisms(self):
        # TODO: This is also the verification
        # TODO: Need column numbers, and make this work for any kind of factor.
        if not hasattr(self, "_morphisms_"):
            return 0
        return 0
