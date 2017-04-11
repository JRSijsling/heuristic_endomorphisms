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

    def calculate_geometric_representations(self):
        if not hasattr(self, "_geo_rep_list_"):
            self._P_ = self.period_matrix()
            self._geo_rep_approx_ = magma.GeometricEndomorphismBasisApproximations(self._P_)
            self._geo_rep_pol_ = magma.RelativeMinimalPolynomialsMatrices(self._geo_rep_approx_[1], self.base_field)
            self._endo_fod_ = Relative_Splitting_Field(self._geo_rep_pol_, bound = self.bound)
            self._geo_rep_list_ = magma.GeometricEndomorphismBasis(self._geo_rep_approx_, self._endo_fod_)
            self._geo_rep_dict_ = dict_rep(self._geo_rep_list_)

    def geometric_representations(self):
        self.calculate_geometric_representations()
        return self._geo_rep_dict_['tangent']

    def endomorphism_field(self):
        self.calculate_geometric_representations()
        return self._endo_fod_

    def geometric(self):
        self.calculate_geometric_representations()
        return OverField(self, K = "geometric")

    def over_base(self):
        self.calculate_geometric_representations()
        return OverField(self, K = "base")

    def over_field(self, K):
        self.calculate_geometric_representations()
        return OverField(self, K = K)

    def lattice(self):
        if not hasattr(self, "_lat_dict_"):
            self.calculate_geometric_representations()
            self._lat_dict_ = dict_lattice(magma.EndomorphismLattice(self._geo_rep_list_))
        return Lattice(self, self._lat_dict_)

    def rosati_involution(self, A):
        self.calculate_geometric_representations()
        return magma.RosatiInvolution(self._geo_rep_list_, A)

    def degree_estimate(self, A):
        self.calculate_geometric_representations()
        return magma.DegreeEstimate(self._geo_rep_list_, A)

    def verify_algebra(self):
        # TODO: Davide
        return True

    def verify_saturated(self):
        self.calculate_geometric_representations()
        return magma.VerifySaturated(self._geo_rep_list_, self._P_)

    def verify_representations(self):
        self.calculate_geometric_representations()
        XL, P0L, AsL = magma.NonWeierstrassBasePointHyp(self.X, self.base_field, self._geo_rep_list_, nvals = 3)
        d = self.degree_estimate(AL)
        return magma.VerifyRepresentations(XL, P0L, AsL, LowerBound = 2*d + 2)

    def verify(self):
        return (self.verify_algebra() and self.verify_saturated() and self.verify_representations())

    def decomposition(self):
        if not hasattr(self, "_lat_dict_"):
            self._lat_dict_ = self.lattice()
        return Decomposition(self)

class Lattice:
    def __init__(self, End, Lat):
        self.X = End.X
        self._lat_dict_ = Lat

    def __repr__(self):
        return repr_lattice(self)

    def representations(self):
        return [ rep['tangent'] for rep in self._lat_dict_['representations'] ]

    def structures(self):
        return self._lat_dict_['algebras']

    def descriptions(self):
        return self._lat_dict_['descriptions']

    def pretty_print(self):
        return pretty_print_lattice(self._lat_dict_['descriptions'], magma.Genus(self.X), 'K', 'x')

class OverField:
    def __init__(self, End, K = "geometric"):
        self.X = End.X
        self.base_field = End.base_field
        self._geo_rep_list_ = End._geo_rep_list_
        self._geo_rep_dict_ = End._geo_rep_dict_
        if K == "geometric":
            self.field = magma.BaseRing(self._geo_rep_dict_['tangent'][1])
        elif K == "base":
            self.field = self.base_field
        else:
            self.field = magma(K)

    def __repr__(self):
        return repr_over_field(self)

    def representations(self):
        if not hasattr(self, "_struct_"):
            self._struct_ = dict_structure(magma.EndomorphismStructure(self._geo_rep_list_, self.field))
        return self._struct_['representation']['tangent']

    # NOTE: Using the Magma functionality EndomorphismAlgebraAndDescription
    # it is straightforward to get versions that do not calculate the Sato-Tate
    # group, which could sometimes give a slight performance gain.
    # Alternatively, it could be a flag whether to include this group or not.
    # Arguably, including Sato-Tate by default violates modularity.
    # Yet for now I see no supremely cogent reason to modify the current
    # approach, which always calculates the Sato-Tate group in one go.
    def algebra(self):
        if not hasattr(self, "_struct_"):
            self._struct_ = dict_structure(magma.EndomorphismStructure(self._geo_rep_list_, self.field))
        return self._struct_['algebra']

    def description(self):
        if not hasattr(self, "_struct_"):
            self._struct_ = dict_structure(magma.EndomorphismStructure(self._geo_rep_list_, self.field))
        return self._struct_['description']

    def pretty_print(self):
        if not hasattr(self, "_struct_"):
            self._struct_ = dict_structure(magma.EndomorphismStructure(self._geo_rep_list_, self.field))
        return pretty_print_over_field(self._struct_['description'], magma.Genus(self.X), 'K')

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
