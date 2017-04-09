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

class EndomorphismData:
    def __init__(self, X, prec, bound = 0, have_oldenburg = False):
        self.X = magma(X)
        self.base_field = magma.BaseRing(self.X)
        self.prec = magma(prec)
        self.bound = magma(bound)
        self.have_oldenburg = magma(have_oldenburg)

    def __repr__(self):
        return repr_endomorphism_data(self)

    def period_matrix(self):
        if not hasattr(self, "_P_"):
            self._eqsCC_ = magma.EmbedCurveEquations(self.X, self.prec)
            self._P_ = magma.PeriodMatrix(self._eqsCC_, HaveOldenburg = self.have_oldenburg)
        return self._P_

    def geometric_representations(self):
        if not hasattr(self, "_geo_reps_"):
            self._P_ = self.period_matrix()
            self._geo_reps_approx_ = magma.GeometricEndomorphismBasisApproximations(self._P_)
            self._geo_reps_pol_ = magma.RelativeMinimalPolynomialsMatrices(self._geo_reps_approx_[1], self.base_field)
            self._endo_fod_ = Relative_Splitting_Field(self._geo_reps_pol_, bound = self.bound)
            self._geo_reps_ = magma.GeometricEndomorphismBasis(self._geo_reps_approx_, self._endo_fod_)
            self._geo_reps_tang_ = self._geo_reps_[1]
            self._geo_reps_hom_ = self._geo_reps_[2]
            self._geo_reps_approx_ = self._geo_reps_[3]
        return self._geo_reps_tang_

    def endomorphism_field(self):
        self._geo_reps_tang_ = self.geometric_representations()
        return self._endo_fod_

    def geometric(self):
        self._geo_reps_tang_ = self.geometric_representations()
        return OverField(self, K = "geometric")

    def over_base(self):
        self._geo_reps_tang_ = self.geometric_representations()
        return OverField(self, K = "base")

    def over_field(self, K):
        self._geo_reps_tang_ = self.geometric_representations()
        return OverField(self, K = K)

    def lattice(self):
        if not hasattr(self, "_lat_"):
            self._geo_reps_tang_ = self.geometric_representations()
            self._lat_ = magma.EndomorphismLattice(self._geo_reps_)
        return Lattice(self, self._lat_)

    def rosati_involution(self, A):
        self._geo_reps_tang_ = self.geometric_representations()
        return magma.RosatiInvolution(self._geo_reps_, A)

    def degree_estimate(self, A):
        self._geo_reps_tang_ = self.geometric_representations()
        return magma.DegreeEstimate(self._geo_reps_, A)

    def verify_algebra(self):
        # TODO: Davide
        return True

    def verify_saturated(self):
        self._geo_reps_tang_ = self.geometric_representations()
        return magma.VerifySaturated(self._geo_reps_, self._P_)

    def verify_representations(self):
        self._geo_reps_tang_ = self.geometric_representations()
        XL, P0L, AsL = magma.NonWeierstrassBasePointHyp(self.X, self.base_field, self._geo_reps_tang_, nvals = 3)
        d = self.degree_estimate(AL)
        return magma.VerifyRepresentations(XL, P0L, AsL, LowerBound = 2*d + 2)

    def verify(self):
        return (self.verify_algebra() and self.verify_saturated() and self.verify_representations())

    def decomposition(self):
        self._P_ = self.period_matrix()
        self._lat_ = self.lattice()
        return Decomposition(self)

class Lattice:
    def __init__(self, End, Lat):
        self.X = End.X
        self._lat_ = Lat
        self._lat_reps_ = Lat[1]
        self._lat_structs_ = Lat[2]
        self._lat_descs_ = Lat[3]

    def __repr__(self):
        return repr_lattice(self)

    def representations(self):
        return [ [ rep[1], rep[2] ] for rep in self._lat_reps_ ]

    def structures(self):
        return self._lat_structs_

    def descriptions(self):
        return self._lat_descs_

    def pretty_print(self):
        return pretty_print_lattice(self._lat_descs_, magma.Genus(self.X), 'F', 'x')

class OverField:
    def __init__(self, End, K = "geometric"):
        self.X = End.X
        self.base_field = End.base_field
        self.field = K
        self._geo_reps_ = End._geo_reps_

    def __repr__(self):
        return repr_over_field(self)

    def representations(self):
        if not hasattr(self, "_reps_"):
            if self.field == "geometric":
                self._reps_ = self._geo_reps_
            elif self.field == "base":
                self._reps_ = magma.EndomorphismBasis(self._geo_reps_, self.base_field)
            else:
                self._reps_ = magma.EndomorphismBasis(self._geo_reps_, magma(self.field))
            self._reps_tang_ = self._reps_[1]
            self._reps_hom_ = self._reps_[2]
            self._reps_approx_ = self._reps_[3]
        return self._reps_tang_

    def structure(self):
        if not hasattr(self, "_struct_"):
            self._reps_tang_ = self.representations()
            self._struct_, self._desc_ = magma.EndomorphismStructure(self._reps_, nvals = 2)
        return self._struct_

    def description(self):
        if not hasattr(self, "_desc_"):
            self._reps_tang_ = self.representations()
            self._struct_, self._desc_ = magma.EndomorphismStructure(self._reps_, nvals = 2)
        return self._desc_

    def pretty_print(self):
        if not hasattr(self, "_desc_"):
            self._reps_tang_ = self.representations()
            self._struct_, self._desc_ = magma.EndomorphismStructure(self._reps_, nvals = 2)
        return pretty_print_over_field(self._desc_, magma.Genus(self.X), 'F')

class Decomposition:
    # We take a smallest field over which everything occurs
    def __init__(self, End):
        self.X = End.X

    def __repr__(self):
        return repr_decomposition(self)

    def field_of_definition(self):
        if not hasattr(self, "_decomp_fod_"):
            return 0
        return 0

    def idempotents(self):
        if not hasattr(self, "_idems_"):
            return 0
        return 0

    def putative_factors(self):
        if not hasattr(self, "_factors_"):
            self._decomp_fod_ = self.field_of_definition()
            self._idems_ = self.idempotents()
            self._factors_ = magma.Factors(self._idems_)
        return self._factors_

    def verified_morphisms(self):
        # TODO: Need column numbers, and make this work for any kind of factor.
        if not hasattr(self, "_morphisms_"):
            return 0
        return 0
