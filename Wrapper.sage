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
        if self.curve_type == "hyperelliptic":
            if self.h == 0:
                return "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(self.f))
            else:
                return "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(self.h), str(self.f))
        elif self.curve_type == "plane":
            return "The endomorphism data of the plane curve over QQ defined by {}".format(str(self.F))

    def period_matrix(self):
        if not hasattr(self, "_P_"):
            if self.curve_type == "hyperelliptic":
                self._P_ = magma.PeriodMatrixHyperelliptic(self._fCC_, self._hCC_, HaveOldenburg = self.have_oldenburg)
            elif self.curve_type == "plane":
                self._P_ = magma.PeriodMatrixPlane(self._FCC_, HaveOldenburg = self.have_oldenburg)
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

    def verify_ring(self):
        self._geo_reps_ = self.geometric_representations()
        if not hasattr(self, "_is_verified_ring_"):
            self._is_verified_ring_ = magma.VerifyRing(self._Rs_, self._P_)
        return self._is_verified_ring_

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

    def geometric_representations_check(self, bound = 2^10):
        self._geo_reps_ = self.geometric_representations()
        XL, P0L, AsL = magma.NonWeierstrassBasePointHyp(self.X, self.base_field, self._AsAlg_, nvals = 3)
        for AL in AsL:
            if magma.IsScalar(AL):
                tests.append(True)
            else:
                # TODO: This transpose should go
                d = self.degree_estimate(AL)
                ALt = magma.Transpose(AL)
                div = magma.CantorMorphismFromMatrixSplit(XL, P0L, ALt, LowerBound = 2*d + 2)
                tests.append(True)
        return all(tests)

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
        # TODO: Small phrase
        return "Some lattice dude!"
        # TODO: Make what follows a conversion function, Galois group
        statement = """Smallest field over which all endomorphisms are defined:\nGalois number field K = QQ (a) with defining polynomial %s\n\n""" % intlist_to_poly(self._frep_)
        for ED in self._lat_descs_:
            statement += """Over subfield F with generator %s with minimal polynomial %s:\n""" % (strlist_to_nfelt(ED[0][1], 'a'), intlist_to_poly(ED[0][0]))
            statement += endo_statement(ED[1], ED[2], ED[3], r'F')
            #statement += st_group_statement(ED[4])
            #statement += gl2_simple_statement(ED[1], ED[2])
            statement += '\n'
        return statement

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
        self.f = EndJac.f
        self.h = EndJac.h
        self.base_field = EndJac.base_field
        self.field = K
        self._geo_reps_ = EndJac._geo_reps_

    def __repr__(self):
        if self.h == 0:
            pre = "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(self.f))
        else:
            pre = "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(self.h), str(self.f))
        if self.field == "geometric":
            post = " over the algebraic closure of its base field"
        elif self.field == "base":
            post = " over its base field"
        else:
            post = " over " + str(self.field)
        return pre + post

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
        # TODO: Add conversion function, a single prettification
        return self._desc_

class Decomposition:
    def __init__(self, EndJac):
        self.f = EndJac.f
        self.h = EndJac.h

    def __repr__(self):
        if self.h == 0:
            return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(self.f))
        else:
            return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(self.h), str(self.f))

    def field_of_definition(self):
        # TODO: Should be calculated with below
        return 0

    def factors(self):
        # TODO: Defer to Magma
        K = self.field_of_definition()
        if not hasattr(self, "_ECs_"):
            Lats = magma.LatticesFromIdempotents(self._idems_[2], self._P_, epscomp = self._epscomp_, epsLLL = self._epsLLL_, epsinv = self._epsinv_)
            self._ECs_rep_ = [ Elliptic_Curve_From_Lattice(Lat.sage(), self._fsubrep_opt_, self._fsubhom_opt_, prec = self.prec, epscomp = self._epscomp_, epsLLL = self._epsLLL_) for Lat in Lats ]
        return [ magma.EllipticCurve([ -K(EC_rep[0])/48, -K(EC_rep[1])/864]).sage() for EC_rep in self._ECs_rep_ ]

    def certificate_g2(self):
        # TODO: Need column numbers, and make this work for any kind of factor. Rest should be like endo verification.
        if len(self._idems_[1]) == 0:
            return [ ]
        return 0
