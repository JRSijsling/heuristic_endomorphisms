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
        self.g = magma.Genus(self.X)
        self.base_field = magma.BaseRing(self.X)
        self.prec = magma(prec)
        self.bound = magma(bound)
        self.have_oldenburg = magma(have_oldenburg)
        self._P_ = self.period_matrix()
        self._calculate_geometric_representations_()

    def __repr__(self):
        return repr_endomorphism_data(self)

    def period_matrix(self):
        self._eqsCC_ = magma.EmbedCurveEquations(self.X, self.prec)
        self._P_ = magma.PeriodMatrix(self._eqsCC_, HaveOldenburg = self.have_oldenburg)
        return self._P_

    def _calculate_geometric_representations_(self):
        self._P_ = self.period_matrix()
        _geo_rep_partial_ = magma.GeometricEndomorphismRepresentationPartial(self._P_)
        _geo_rep_pol_ = magma.RelativeMinimalPolynomialsPartial(_geo_rep_partial_, self.base_field)
        self._endo_fod_ = Relative_Splitting_Field(_geo_rep_pol_, bound = self.bound)
        self._geo_rep_list_ = magma.GeometricEndomorphismRepresentationRecognition(_geo_rep_partial_, self._endo_fod_)
        self._geo_rep_dict_ = dict_rep(self._geo_rep_list_)

    def endomorphism_field(self):
        return self._endo_fod_

    def geometric(self):
        return OverField(self, K = "geometric")

    def over_base(self):
        return OverField(self, K = "base")

    def over_field(self, K):
        return OverField(self, K = K)

    def lattice(self):
        if not hasattr(self, "_lat_"):
            self._lat_ = Lattice(self)
        return self._lat_

    def decomposition(self):
        if not hasattr(self, "_lat_"):
            self._lat_ = self.lattice()
        return Decomposition(self)

    def rosati_involution(self, A):
        return magma.RosatiInvolution(self._geo_rep_list_, A)

    def degree_estimate(self, A):
        return magma.DegreeEstimate(self._geo_rep_list_, A)

    def dimension_algebra(self):
        return len(self._geo_rep_list_)

    def verify_algebra(self):
        # TODO: Integrate over Davide and Edgar
        self._test_alg_ = True
        return self._test_alg_

    def verify_saturated(self):
        self._sat_test_, self._sat_cert_ =  magma.VerifySaturated(self._geo_rep_list_, self._P_, nvals = 2)
        return self._sat_test_

    def set_base_point(self):
        if not hasattr(self, "base_point"):
            self.base_point = magma.NonWeierstrassBasePoint(self.X, self._endo_fod_)

    def correspondence(self, A):
        # TODO: Add bounds
        self.set_base_point()
        test, cert = magma.Correspondence(self.X, self.base_point, self.X, self.base_point, A, nvals = 2)
        return cert

    def verify_representations(self):
        self._rep_test_ = True
        for gen in self._geo_rep_dict_:
            genTan = gen['tangent']
            corresp = self.correspondence(genTan)
            if not corresp:
                self._rep_test_ = False
            else:
                gen['corresp'] = corresp
        return self._rep_test_

    def verify(self):
        return (self.verify_algebra() and self.verify_saturated() and self.verify_representations())

class OverField:
    def __init__(self, Endo, K = "geometric"):
        self.X = Endo.X
        self.g = Endo.g
        self.base_field = Endo.base_field
        self._P_ = Endo._P_
        self._geo_rep_list_ = Endo._geo_rep_list_
        if K == "geometric":
            self.field = Endo.endomorphism_field()
        elif K == "base":
            self.field = self.base_field
        else:
            self.field = magma(K)
        self._list_ = magma.EndomorphismStructure(self._geo_rep_list_, self.field)
        self._desc_ = desc_structure(self._list_)

    def __repr__(self):
        return repr_over_field(self)

    def pretty_print(self):
        return pretty_print_over_field_description(self._desc_, self.g, 'K')

    def _calculate_dictionary_(self):
        if not hasattr(self, "_dict_"):
            self._dict_ = dict_structure(self._list_)

    def full(self):
        self._calculate_dictionary_()
        return self._dict_

    def representation(self):
        self._calculate_dictionary_()
        return magma([ gen['tangent'] for gen in self._dict_['representation'] ])

    def algebra(self):
        self._calculate_dictionary_()
        return self._dict_['algebra']

    def description(self):
        self._calculate_dictionary_()
        return self._dict_['description']

    def rosati_involution(self, A):
        return magma.RosatiInvolution(self._list_[_index_dict_['representation']], A)

    def degree_estimate(self, A):
        return magma.DegreeEstimate(self._list_[_index_dict_['representation']], A)

    def dimension_algebra(self):
        return len(self._list_)

    def verify_algebra(self):
        # TODO: Integrate over Davide and Edgar
        self._test_alg_ = True
        return self._test_alg_

    def verify_saturated(self):
        self._sat_test_, self._sat_cert_ =  magma.VerifySaturated(self._list_[1], self._P_, nvals = 2)
        return self._sat_test_

    def set_base_point(self):
        if not hasattr(self, "base_point"):
            self.base_point = magma.NonWeierstrassBasePoint(self.X, self.field)

    def correspondence(self, A):
        self.set_base_point()
        test, cert = magma.Correspondence(self.X, self.base_point, self.X, self.base_point, A, nvals = 2)
        return cert

    def verify_representations(self):
        self._calculate_dictionary_()
        self._rep_test_ = True
        for gen in self._dict_['representation']:
            genTan = gen['tangent']
            corresp = self.correspondence(genTan)
            if not corresp:
                self._rep_test_ = False
            else:
                gen['corresp'] = corresp
        return self._rep_test_

    def verify(self):
        return (self.verify_algebra() and self.verify_saturated() and self.verify_representations())

class Lattice:
    def __init__(self, Endo):
        self.X = Endo.X
        self.g = Endo.g
        self._geo_rep_list_ = Endo._geo_rep_list_
        self._list_ = magma.EndomorphismLattice(self._geo_rep_list_)
        self._desc_ = desc_lattice(self._list_)

    def __repr__(self):
        return repr_lattice(self)

    def _calculate_dictionary_(self):
        if not hasattr(self, "_dict_"):
            self._dict_ = dict_lattice(self._list_)

    def full(self):
        self._calculate_dictionary_()
        return self._dict_

    def representations(self):
        self._calculate_dictionary_()
        list_to_fill = [ ]
        for dict_pair in self._dict_:
            dict_to_fill = dict()
            dict_to_fill['field'] = dict_pair['field']['magma']
            dict_to_fill['representation'] =  magma([ gen['tangent'] for gen in dict_pair['structure']['representation'] ])
            list_to_fill.append(dict_to_fill)
        return list_to_fill

    def algebras(self):
        self._calculate_dictionary_()
        list_to_fill = [ ]
        for dict_pair in self._dict_:
            dict_to_fill = dict()
            dict_to_fill['field'] = dict_pair['field']['magma']
            dict_to_fill['algebra'] = dict_pair['structure']['algebra']
            list_to_fill.append(dict_to_fill)
        return list_to_fill

    def descriptions(self):
        self._calculate_dictionary_()
        list_to_fill = [ ]
        for dict_pair in self._dict_:
            dict_to_fill = dict()
            dict_to_fill['field'] = dict_pair['field']['magma']
            dict_to_fill['description'] = dict_pair['structure']['description']
            list_to_fill.append(dict_to_fill)
        return list_to_fill

    def pretty_print(self):
        return pretty_print_lattice_description(self._desc_, self.g, 'K', 'x')

class Decomposition:
    def __init__(self, Endo):
        self.X = Endo.X
        self.g = Endo.g
        self._P_ = Endo._P_
        idems, self.field = magma.IdempotentsFromLattice(Endo._lat_._list_, nvals = 2)
        self._facs_ = [ ]
        for idem in idems:
            fac = dict()
            fac['field'] = self.field
            fac['idem'] = dict_gen(idem)
            lat, proj = magma.ProjectionFromIdempotent(self._P_, idem, nvals = 2)
            fac['proj'] = dict_gen(proj);
            fac['factor'] = { 'analytic': lat }
            self._facs_.append(fac)

    def __repr__(self):
        return repr_decomposition(self)

    def full(self):
        return self._facs_

    def idempotents(self):
        return [ fac['idem']['tangent'] for fac in self._facs_ ]

    def projections(self):
        return [ fac['proj']['tangent'] for fac in self._facs_ ]

    def _calculate_factors_(self):
        for fac in self._facs_:
            if not 'algebraic' in fac['factor'].keys():
                fac['factor']['algebraic'] = magma.FactorReconstruct(fac['factor']['analytic'], fac['field'])

    def factors(self):
        self._calculate_factors_()
        return [ fac['factor']['algebraic'] for fac in self._facs_ ]

    def set_base_point(self):
        if not hasattr(self, "base_point"):
            self.base_point = magma.NonWeierstrassBasePoint(self.X, self.field)

    def set_base_point_factor(self, fac):
        if not 'base_point' in fac['factor'].keys():
            Y = fac['factor']['algebraic']
            K = self.field
            fac['factor']['base_point'] = magma.NonWeierstrassBasePoint(Y, K)

    def correspondence(self, fac):
        self.set_base_point()
        self.set_base_point_factor(fac)
        P = self.base_point
        Y = fac['factor']['algebraic']
        Q = fac['factor']['base_point']
        A = fac['proj']['tangent']
        test, cert = magma.Correspondence(self.X, P, Y, Q, A, nvals = 2)
        return cert

    def verify(self):
        self._facs_test_ = True
        for fac in self._facs_:
            corresp = self.correspondence(fac)
            if not corresp:
                self._facs_test_ = False
            else:
                fac['proj']['corresp'] = corresp
        return self._facs_test_
