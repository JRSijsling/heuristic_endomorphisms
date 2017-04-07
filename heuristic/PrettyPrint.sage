# TODO: Avoid boilerplate by ReprCurve

def ReprEndomorphismData(End):
    if End.curve_type == "hyperelliptic":
        if End.h == 0:
            return "The endomorphism data of the hyperelliptic curve y^2 = {}".format(str(End.f))
        else:
            return "The endomorphism data of the hyperelliptic curve y^2 + ({})*y = {}".format(str(End.h), str(End.f))
    elif End.curve_type == "plane":
        return "The endomorphism data of the plane curve {} = 0".format(str(End.F))

def ReprLattice(Lat):
    # TODO: Small phrase
    return "Some lattice dude!"
    # TODO: Make what follows a conversion function, Galois group
    statement = """Smallest field over which all endomorphisms are defined:\nGalois number field K = QQ (a) with defining polynomial %s\n\n""" % intlist_to_poly(Lat._frep_)
    for ED in Lat._lat_descs_:
        statement += """Over subfield F with generator %s with minimal polynomial %s:\n""" % (strlist_to_nfelt(ED[0][1], 'a'), intlist_to_poly(ED[0][0]))
        statement += endo_statement(ED[1], ED[2], ED[3], r'F')
        #statement += st_group_statement(ED[4])
        #statement += gl2_simple_statement(ED[1], ED[2])
        statement += '\n'
    return statement

def ReprOverField(OF):
    if OF.h == 0:
        pre = "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(OF.f))
    else:
        pre = "The endomorphism data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(OF.h), str(OF.f))
    if OF.field == "geometric":
        post = " over the algebraic closure of its base field"
    elif OF.field == "base":
        post = " over its base field"
    else:
        post = " over " + str(OF.field)
    return pre + post

def ReprDescription(desc):
    return "boink"

def ReprDecomposition(decomp):
    if decomp.h == 0:
        return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 = {}".format(str(decomp.f))
    else:
        return "The decomposition data of the hyperelliptic curve over QQ defined by y^2 + ({})*y = {}".format(str(decomp.h), str(decomp.f))
