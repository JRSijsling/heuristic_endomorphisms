def CurveType(X):
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

def ReprCurve(X):
    curve_type = CurveType(X)
    if curve_type == "hyperelliptic":
        f, h = X.hyperelliptic_polynomials()
        if h == 0:
            return " the hyperelliptic curve y^2 = {}".format(str(f))
        else:
            return " the hyperelliptic curve y^2 + ({})*y = {}".format(str(h), str(f))
    elif End.curve_type == "plane":
        F = X.defining_polynomial()
        return " the plane curve {} = 0".format(str(F))

def ReprEndomorphismData(End):
    return "The endomorphism data of" + ReprCurve(End.X)

def ReprLattice(Lat):
    return "The endomorphism lattice of" + ReprCurve(Lat.X)
    # TODO: Make what follows a conversion function, Galois group
    statement = """Smallest field over which all endomorphisms are defined:\nGalois number field K = QQ (a) with defining polynomial %s\n\n""" % intlist_to_poly(Lat._frep_)
    for ED in Lat._lat_descs_:
        statement += """Over subfield F with generator %s with minimal polynomial %s:\n""" % (strlist_to_nfelt(ED[0][1], 'a'), intlist_to_poly(ED[0][0]))
        statement += endo_statement(ED[1], ED[2], ED[3], r'F')
        #statement += st_group_statement(ED[4])
        #statement += gl2_simple_statement(ED[1], ED[2])
        statement += '\n'
    return statement

def ReprOverField(overfield):
    pre = "The endomorphism structure of" + ReprCurve(overfield.X)
    if overfield.field == "geometric":
        post = " over the algebraic closure of its base field"
    elif overfield.field == "base":
        post = " over its base field"
    else:
        post = " over " + str(overfield.field)
    return pre + post

def ReprDescription(desc):
    return "boink"

def ReprDecomposition(decomp):
    return "The decomposition structure of" + ReprCurve(decomp.X)
