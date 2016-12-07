"""
 *  Class wrappers for the algorithms
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
"""

class OverField:
    def __init__(self, EndJac, K = "geometric"):
        self.g = EndJac.g
        self.h = EndJac.h
        self._geo_reps_ = EndJac._geo_reps_
        self._lat_ = EndJac._lat_
        self.field = K

    def __repr__(self):
        if self.h == 0:
            pre = "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(self.g))
        else:
            pre = "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(self.h), str(self.g))
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

class Lattice:
    def __init__(self, EndJac):
        self.g = EndJac.g
        self.h = EndJac.h
        self._frep_ = EndJac._frep_
        self._lat_ = EndJac._lat_

    def __repr__(self):
        statement = ''
        if self.h == 0:
            statement += "The endomorphism lattice of the hyperelliptic curve over QQ defined by y^2 = {}\n\n".format(str(self.g))
        else:
            statement += "The endomorphism lattice of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}\n\n".format(str(self.h), str(self.g))
        statement += """Smallest field over which all endomorphisms are defined:\nGalois number field K = QQ (a) with defining polynomial %s\n\n""" % intlist_to_poly(self._frep_)
        for ED in self._lat_:
            statement += """Over subfield F with generator %s with minimal polynomial %s:\n""" % (strlist_to_nfelt(ED[0][1], 'a'), intlist_to_poly(ED[0][0]))
            statement += endo_statement(ED[1], ED[2], ED[3], r'F')
            statement += st_group_statement(ED[4])
            statement += gl2_simple_statement(ED[1], ED[2])
            statement += '\n'
        return statement

class Decomposition:
    def __init__(self, EndJac):
        self.g = EndJac.g
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
            return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(self.g))
        else:
            return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(self.h), str(self.g))

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
        spl_fod_data = Canonize_Subfield_And_Idempotents(self._fsubgen_.sage(), self._frep_, self._idems_)
        self._spl_fod_ = spl_fod_data
        self._fsubrep_opt_ = spl_fod_data[0]
        self._fsubgen_opt_ = spl_fod_data[1]
        self._fsubhom_opt_ = magma.InducedEmbedding(self._fsubgen_opt_, self._fhom_, self._frep_).sage()
        K = self.field_of_definition()
        Lats, col_numbers = magma.LatticesFromIdempotents(self._idems_[2], self._P_, epscomp = self._epscomp_, epsLLL = self._epsLLL_, epsinv = self._epsinv_, nvals = 2)
        self._ECs_rep_ = [ Elliptic_Curve_From_Lattice(Lat.sage(), self._fsubrep_opt_, self._fsubhom_opt_, prec = self.prec, epscomp = self._epscomp_, epsLLL = self._epsLLL_) for Lat in Lats ]
        # Convert everything and pass to Magma
        fX = magma(4*self.g + self.h^2)
        X = magma.HyperellipticCurve(fX)
        As = magma(spl_fod_data[2])
        As = magma.ProjectFromColumnNumbers(As, col_numbers)
        XL, P0L, AsL, L = magma.NonWeierstrassBasePointHyp(X, K, As, nvals = 4)
        EsL, Q0sL = magma.EllipticCurvesFromRepresentations(magma(self._ECs_rep_), K, L, nvals = 2)
        # TODO: Next line not used yet, should be upper bound
        degs = magma([ magma.DecompositionDegree(idem) for idem in self._idems_[3] ])
        zipped_list = zip(EsL, Q0sL, AsL)
        return [ magma.CantorMorphismFromMatrixSplit(XL, P0L, entry[0], entry[1], entry[2]) for entry in zipped_list ]

class EndomorphismData:
    def __init__(self, g, h = 0, prec = prec):
        RQQ.<x> = PolynomialRing(QQ)
        self.g = RQQ(g)
        self.h = RQQ(h)
        self.prec = prec
        self._epscomp_ = 10^(-self.prec + 30)
        self._epsLLL_ = 5^(-self.prec + 7)
        self._epsinv_ = 2^(-10)
        self._Bound_ = Bound
        self._lat_ = None

    def __repr__(self):
        if self.h == 0:
            return "The functor of endomorphism data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(self.g))
        else:
            return "The functor of endomorphism data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(self.h), str(self.g))

    def period_matrix(self):
        if not hasattr(self, "_P_"):
            self._P_ = magma.PeriodMatrix(self.h, self.g, prec = self.prec)
        return self._P_.sage()

    def field_of_definition(self):
        if not hasattr(self, "_frep_"):
            geo_reps = self.geometric_representations()
        RQQ.<x> = PolynomialRing(QQ)
        K.<r> = NumberField(RQQ(self._frep_))
        return K

    def geometric_representations(self):
        if not hasattr(self, "_geo_reps_"):
            P = self.period_matrix()
            self._geo_reps_app_ = magma.GeometricEndomorphismBasisFromPeriodMatrix(self._P_, epscomp = self._epscomp_, epsLLL = self._epsLLL_, epsinv = self._epsinv_, nvals = 2)
            self._geo_reps_pol_ = [ magma.PolynomializeMatrix(A, epscomp = self._epscomp_, epsLLL = self._epsLLL_) for A in self._geo_reps_app_[0] ]
            self._frep_ = Common_Splitting_Field(self._geo_reps_pol_, Bound = self._Bound_)
            L = magma.AlgebraizeMatricesInField(self._geo_reps_app_[0], self._geo_reps_pol_, self._frep_, epscomp = self._epscomp_, nvals = 2)
            self._geo_reps_alg_ = L[0]
            self._fhom_ = L[1]
            self._geo_reps_ = [ self._geo_reps_alg_, self._geo_reps_app_[0], self._geo_reps_app_[1] ]
        # For some weird reason this fails:
        #return [ [ rep.sage() for rep in reps ] for reps in self._geo_reps_ ]
        return self._geo_reps_

    def geometric_representations_check(self, bound = 2^10):
        # TODO: Split case
        X = magma.HyperellipticCurve(4*self.g + self.h^2)
        As = self.geometric_representations()[0]
        K = magma.BaseRing(magma.Parent(As[1]))
        X, P0, AsL = magma.NonWeierstrassBasePointHyp(X, K, As, nvals = 3)
        tests = [ ]
        for i in [1..len(AsL)]:
            if magma.IsScalar(As[i]):
                tests.append(True)
            else:
                d = self.degree_estimate(As[i])
                tests.append(d)
                AtL = magma.Transpose(AsL[i])
                div = magma.CantorMorphismFromMatrixSplit(X, P0, AtL, DegreeBound = 2*d + 2)
                # TODO: This justs appends True for now because in fact the current test will never stop if there is an error... TBD
                tests.append(True)
        return all(tests)

    def geometric(self):
        geo_reps = self.geometric_representations()
        return OverField(self, K = "geometric")

    def over_base(self):
        geo_reps = self.geometric_representations()
        return OverField(self, K = "base")

    def over_field(self, K):
        geo_reps = self.geometric_representations()
        return OverField(self, K = K)

    def rosati_involution(self, A):
        geo_reps = self.geometric_representations()
        return magma.RosatiInvolution(geo_reps[0], geo_reps[1], geo_reps[2], A)

    def degree_estimate(self, A):
        geo_reps = self.geometric_representations()
        return magma.DegreeEstimate(geo_reps[0], geo_reps[1], geo_reps[2], A)

    def lattice(self):
        if not self._lat_:
            AsAlg, As, Rs = self.geometric_representations()
            L = magma.EndomorphismLatticeG2(AsAlg, As, Rs, Geometric = false, AddTensor = true, AddRing = true, AddSatoTate = true, AddDecomposition = true, nvals = 3)
            self._lat_ = EDs_sagified(L[0], self._frep_)
            self._fsubgen_ = L[1]
            self._idems_ = L[2]
        return Lattice(self)

    def decomposition(self):
        P = self.period_matrix()
        lat = self.lattice()
        return Decomposition(self)
