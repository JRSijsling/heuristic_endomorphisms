"""
 *  Class wrappers for the algorithms
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
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
            self._fod_ = magma.BaseRing(self.f)
            embedded_list = magma.EmbedAsComplexPolynomials([self.f, self.h], prec)
            self._fCC_, self._hCC_ = embedded_list
        elif self.curve_type == "plane":
            self.F = magma(X.defining_polynomial())
            self._fod_ = magma.BaseRing(self.F)
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
            self._AsPol_ = magma.RelativeMinimalPolynomialsMatrices(self._geo_reps_[1], self._fod_)
            self._endo_fod_ = Relative_Splitting_Field(self._AsPol_, bound = self.bound)
            return self._endo_fod_
            self._geo_reps_ = magma.GeometricEndomorphismBasisRepresentations(self._geo_reps_, self._fod_)
        return self._geo_reps_

    def base_field(self):
        self._geo_reps_ = self.geometric_representations()
        return self._fod_

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

    def lattice_description(self):
        if not hasattr(self, "_lat_"):
            self._geo_reps_ = self.geometric_representations()
            self._lat_ = magma.EndomorphismLatticeDescription(self._geo_reps_)
        return Lattice(self)

#    def rosati_involution(self, A):
#        geo_reps = self.geometric_representations()
#        return magma.RosatiInvolution(geo_reps[0], geo_reps[1], geo_reps[2], A)
#
#    def degree_estimate(self, A):
#        geo_reps = self.geometric_representations()
#        return magma.DegreeEstimate(geo_reps[0], geo_reps[1], geo_reps[2], A)

#    def geometric_representations_check(self, bound = 2^10):
#        # TODO: Split case
#        X = magma.HyperellipticCurve(4*self.f + self.h^2)
#        As = self.geometric_representations()[0]
#        K = magma.BaseRing(magma.Parent(As[1]))
#        X, P0, AsL = magma.NonWeierstrassBasePointHyp(X, K, As, nvals = 3)
#        tests = [ ]
#        for i in [1..len(AsL)]:
#            if magma.IsScalar(As[i]):
#                tests.append(True)
#            else:
#                d = self.degree_estimate(As[i])
#                tests.append(d)
#                AtL = magma.Transpose(AsL[i])
#                div = magma.CantorMorphismFromMatrixSplit(X, P0, AtL, LowerBound = 2*d + 2)
#                # TODO: This justs appends True for now because in fact the current test will never stop if there is an error... TBD
#                tests.append(True)
#        return all(tests)

#    def decomposition(self):
#        P = self.period_matrix()
#        lat = self.lattice()
#        return Decomposition(self)

class Lattice:
    def __init__(self, EndJac):
        # TODO: Make this different
        self._frep_ = EndJac._frep_
        self._lat_ = EndJac._lat_

    def __repr__(self):
        # TODO: Better represetation postponed
        return str(self._lat_)
        statement = """Smallest field over which all endomorphisms are defined:\nGalois number field K = QQ (a) with defining polynomial %s\n\n""" % intlist_to_poly(self._frep_)
        for ED in self._lat_:
            statement += """Over subfield F with generator %s with minimal polynomial %s:\n""" % (strlist_to_nfelt(ED[0][1], 'a'), intlist_to_poly(ED[0][0]))
            statement += endo_statement(ED[1], ED[2], ED[3], r'F')
            #statement += st_group_statement(ED[4])
            #statement += gl2_simple_statement(ED[1], ED[2])
            statement += '\n'
        return statement

class OverField:
    def __init__(self, EndJac, K = "geometric"):
        self.f = EndJac.f
        self.h = EndJac.h
        self._geo_reps_ = EndJac._geo_reps_
        self._lat_ = EndJac._lat_
        self.field = K

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
                AsAlg, As, Rs = self._geo_reps_
                self._reps_ = magma.EndomorphismBasisOverField(magma.Rationals(), AsAlg, As, Rs)
            else:
                AsAlg, As, Rs = self._geo_reps_
                self._reps_ = magma.EndomorphismBasisOverField(magma(self.field), AsAlg, As, Rs)
        return self._reps_

    def description(self):
        # TODO: Should be deferred elsewhere
        return 0
        if not hasattr(self, "_desc_"):
            if self._lat_ != None:
                if self.field == "geometric":
                    ED = self._lat_[-1]
                elif self.field == "base":
                    ED = self._lat_[0]
                else:
                    for EDtry in self._lat_[::-1]:
                        L = magma.NumberField(magma.Polynomial(EDtry[0][0]))
                        if magma.IsSubfield(L, magma(K)):
                            ED = EDtry
                            break
                self._desc_ = endo_statement(ED[1], ED[2], ED[3], r'F') + st_group_statement(ED[4]) + gl2_simple_statement(ED[1], ED[2])
            else:
                AsAlg, As, Rs = self._geo_reps_
                if self.field == "geometric":
                    ED = magma.EndomorphismLatticeG2SingleElement(magma.BaseRing(AsAlg[1]), AsAlg, As, Rs)
                elif self.field == "base":
                    ED = magma.EndomorphismLatticeG2SingleElement(magma.Rationals(), AsAlg, As, Rs)
                else:
                    ED = magma.EndomorphismLatticeG2SingleElement(magma(K), AsAlg, As, Rs)
                ED = ED_sagified(ED)
                self._desc_ = endo_statement(ED[0], ED[1], ED[2], r'F') + st_group_statement(ED[3]) + gl2_simple_statement(ED[0], ED[1])
        return self._desc_

class Decomposition:
    def __init__(self, EndJac):
        self.f = EndJac.f
        self.h = EndJac.h
        self._frep_ = EndJac._frep_
        self._fhom_ = EndJac._fhom_
        self._fsubgen_ = EndJac._fsubgen_
        self._idems_ = EndJac._idems_
        self._P_ = EndJac._P_
        self.prec = EndJac.prec
        self._epscomp_ = EndJac._epscomp_
        self._epsLLL_ = EndJac._epsLLL_
        self._epsinv_ = EndJac._epsinv_

    def __repr__(self):
        if self.h == 0:
            return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(self.f))
        else:
            return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(self.h), str(self.f))

    def field_of_definition(self):
        if not hasattr(self, "_spl_fod_"):
            spl_fod_data = Canonize_Subfield(self._fsubgen_.sage(), self._frep_)
            self._spl_fod_ = spl_fod_data
            self._fsubrep_opt_ = spl_fod_data[0]
            self._fsubgen_opt_ = spl_fod_data[1]
            self._fsubhom_opt_ = magma.InducedEmbedding(self._fsubgen_opt_, self._fhom_, self._frep_).sage()
        return magma.NumberField(magma.Polynomial(self._fsubrep_opt_))

    def factors(self):
        K = self.field_of_definition()
        if not hasattr(self, "_ECs_"):
            Lats = magma.LatticesFromIdempotents(self._idems_[2], self._P_, epscomp = self._epscomp_, epsLLL = self._epsLLL_, epsinv = self._epsinv_)
            self._ECs_rep_ = [ Elliptic_Curve_From_Lattice(Lat.sage(), self._fsubrep_opt_, self._fsubhom_opt_, prec = self.prec, epscomp = self._epscomp_, epsLLL = self._epsLLL_) for Lat in Lats ]
        return [ magma.EllipticCurve([ -K(EC_rep[0])/48, -K(EC_rep[1])/864]).sage() for EC_rep in self._ECs_rep_ ]

    def certificate_g2(self):
        if len(self._idems_[1]) == 0:
            return [ ]
        spl_fod_data = Canonize_Subfield_And_Idempotents(self._fsubgen_.sage(), self._frep_, self._idems_)
        self._spl_fod_ = spl_fod_data
        self._fsubrep_opt_ = spl_fod_data[0]
        self._fsubgen_opt_ = spl_fod_data[1]
        self._fsubhom_opt_ = magma.InducedEmbedding(self._fsubgen_opt_, self._fhom_, self._frep_).sage()
        Lats, col_numbers = magma.LatticesFromIdempotents(self._idems_[2], self._P_, epscomp = self._epscomp_, epsLLL = self._epsLLL_, epsinv = self._epsinv_, nvals = 2)
        self._ECs_rep_ = [ Elliptic_Curve_From_Lattice(Lat.sage(), self._fsubrep_opt_, self._fsubhom_opt_, prec = self.prec, epscomp = self._epscomp_, epsLLL = self._epsLLL_) for Lat in Lats ]
        fX = magma(4*self.f + self.h^2)
        X = magma.HyperellipticCurve(fX)
        As = magma.ProjectFromColumnNumbers(spl_fod_data[2], col_numbers)
        K = magma.BaseRing(As[1])
        XL, P0L, AsL, L = magma.NonWeierstrassBasePointHyp(X, K, As, nvals = 4)
        EsL, Q0sL = magma.EllipticCurvesFromRepresentations(magma(self._ECs_rep_), K, L, nvals = 2)
        # TODO: Next line not used yet, should be upper bound
        degs = magma([ magma.DecompositionDegree(idem) for idem in self._idems_[3] ])
        zipped_list = zip(EsL, Q0sL, AsL)
        return [ magma.CantorMorphismFromMatrixSplit(XL, P0L, entry[0], entry[1], entry[2]) for entry in zipped_list ]
