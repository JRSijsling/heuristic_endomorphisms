"""
 *  Representation and pretty print functions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

def repr_curve(X):
    curve_type = magma.CurveType(X)
    if str(curve_type) == "hyperelliptic":
        f, h = magma.HyperellipticPolynomials(X, nvals = 2)
        if magma.IsZero(h):
            return " the hyperelliptic curve y^2 = {}".format(str(f))
        else:
            return " the hyperelliptic curve y^2 + ({})*y = {}".format(str(h), str(f))
    elif str(curve_type) == "plane":
        F = magma.DefiningPolynomial(X)
        return " the plane curve {} = 0".format(str(F))

def repr_endomorphism_data(End):
    return "The endomorphism data of" + repr_curve(End.X)

def repr_lattice(Lat):
    return "The endomorphism lattice of" + repr_curve(Lat.X)

def repr_over_field(over_field):
    pre = "The endomorphism structure of" + repr_curve(over_field.X)
    if over_field.field == "geometric":
        post = " over the algebraic closure of its base field"
    elif over_field.field == "base":
        post = " over its base field"
    else:
        post = " over " + str(over_field.field)
    return pre + post

def repr_decomposition(decomp):
    return "The decomposition structure of" + repr_curve(decomp.X)

def pretty_print_over_field(desc, genus, fieldstring):
    desc_sage = sagify_description(desc)
    desc_dict = dict_description(desc_sage)
    return statements_all(desc_dict, genus, fieldstring)

def pretty_print_lattice(desc_lat, genus, fieldstring, varstring):
    desc_lat = sagify_description(desc_lat)
    statements = [ ]
    for desc in desc_lat:
        statement = ''
        field_desc = desc.pop(0)
        statement += "Over subfield %s:\n" % pretty_print_field(field_desc)
        desc_dict = dict_description(desc)
        statement += statements_all(desc_dict, genus, fieldstring)
        statements.append(statement)
    return '\n\n'.join(statements)

def pretty_print_polynomial_list(s, varstring):
    return str(PolynomialRing(QQ, varstring)(s))

def pretty_print_ring(field_desc, index):
    # TODO: General base field, and for example cyclotomic extensions
    if len(field_desc) == 2:
        return 'ZZ'
    elif len(field_desc) == 3:
        c ,b, a = field_desc
        disc = b**2 - 4*a*c
        if index == 1:
            if disc % 4 == 0:
                return 'ZZ [sqrt(%s)]'% str(disc//4)
            return 'ZZ [(1 + sqrt(%s))/2]' % str(disc)
        if disc % 4 == 0:
            return 'ZZ [%s sqrt(%s)]' % (str(index), str(disc//4))
        if index % 2 == 0:
            return 'ZZ [%s sqrt(%s)]' % (str(index//2), str(disc))
        return 'ZZ [(1 + sqrt(%s))/2]' % (str(index), str(disc))
    else:
        if index == 1:
            return 'Int (%s)' % pretty_print_field(field_desc)
        else:
            return 'Sub (%s, %s)' % (pretty_print_field(field_desc), index)

def pretty_print_field(field_desc):
    # Only makes a difference for at most quadratic fields
    # TODO: General base field, and for example cyclotomic extensions
    if len(field_desc) == 2:
        return 'QQ'
    if len(field_desc) == 3:
        c, b, a = field_desc
        D = b**2 - 4*a*c
        return 'QQ (sqrt(%s))' % D.squarefree_part()
    else:
        return 'QQ [x] / (%s)' % pretty_print_polynomial_list(field_desc)

def sagify_description(desc):
    if str(magma.Type(desc)) == "RngIntElt":
        return ZZ(desc)
    elif str(magma.Type(desc)) == "MonStgElt":
        return str(desc)
    else:
        return [ sagify_description(desc[i]) for i in [1..len(desc)] ]

def dict_factor_QQ(factor_QQ):
    dict_to_fill = dict()
    dict_to_fill['albert_type'] = factor_QQ[0]
    dict_to_fill['base_field'] = factor_QQ[1]
    dict_to_fill['dim_sqrt'] = factor_QQ[2]
    dict_to_fill['disc'] = factor_QQ[3]
    return dict_to_fill

def dict_desc_ZZ(desc_ZZ):
    dict_to_fill = dict()
    dict_to_fill['index'] = desc_ZZ[0]
    dict_to_fill['is_eichler'] = desc_ZZ[1]
    return dict_to_fill

def dict_desc_RR(desc_RR):
    dict_to_fill = dict()
    dict_to_fill['factors'] = desc_RR
    return dict_to_fill

def dict_description(desc):
    dict_to_fill = dict()
    dict_to_fill['factors_QQ'] = [ dict_factor_QQ(factor_QQ) for factor_QQ in desc[0] ]
    dict_to_fill['desc_ZZ'] = dict_desc_ZZ(desc[1])
    dict_to_fill['desc_RR'] = dict_desc_RR(desc[2])
    dict_to_fill['sato_tate'] = ""
    # TODO: Deal with next line and the substitutions that it requires.
    #dict_to_fill['sato_tate'] = desc[3]
    return dict_to_fill

def statements_all(desc_dict, genus, fieldstring):
    statements= [ ]
    statements.append(statement_endomorphisms_QQ(desc_dict, genus, fieldstring) + ' ' + statement_cm(desc_dict, genus, fieldstring))
    statements.append(statement_endomorphisms_ZZ(desc_dict, genus, fieldstring) + ' ' + statement_eichler(desc_dict, genus, fieldstring))
    statements.append(statement_endomorphisms_RR(desc_dict, genus, fieldstring))
    statements.append(statement_sato_tate_group(desc_dict, genus, fieldstring))
    statements.append(statement_gl2(desc_dict, genus, fieldstring) + '; ' + statement_simple(desc_dict, genus, fieldstring))
    return '\n'.join(statements)

def statement_endomorphisms_QQ(desc_dict, genus, fieldstring):
    factor_dicts = desc_dict['factors_QQ']
    statements = [ statement_factor_QQ(factor_dict) for factor_dict in factor_dicts ]
    statement =  ' x '.join(statements)
    return "End (J_%s) ox QQ: " % fieldstring + statement

def statement_endomorphisms_ZZ(desc_dict, genus, fieldstring):
    factors_QQ = desc_dict['factors_QQ']
    desc_ZZ = desc_dict['desc_ZZ']
    if desc_ZZ['index'] == 1:
        statements = [ statement_factor_ZZ_maximal(factor_QQ, desc_ZZ, fieldstring) for factor_QQ in factors_QQ ]
        statement = ' x '.join(statements)
    else:
        statement = statement_factors_ZZ_index(factors_QQ, desc_ZZ, fieldstring)
    return "End (J_%s): " % fieldstring + statement

def statement_endomorphisms_RR(desc_dict, genus, fieldstring):
    return "End (J_%s) ox RR: %s" % (fieldstring, ' x '.join(desc_dict['desc_RR']['factors']))

def statement_factor_QQ(factor_QQ):
    base_str = pretty_print_field(factor_QQ['base_field'])
    d = factor_QQ['dim_sqrt']
    d_str = str(d)
    disc = factor_QQ['disc']
    disc_str = str(disc)

    if factor_QQ['albert_type'] == 'I':
        statement = base_str

    elif factor_QQ['albert_type'] == 'II':
        if disc == 1:
            statement = "M_%s (%s)" % (d_str, base_str)
        elif d == 2:
            statement = "IndefQuat (%s, %s)"  % (base_str, disc_str)
        else:
            statement = "IndefAlg_%s (%s, %s)"  % (d_str, base_str, disc_str)

    elif factor_QQ['albert_type'] == 'III':
        if d == 2:
            statement = "DefQuat (%s, %s)"  % (base_str, disc_str)
        else:
            statement = "DefAlg_%s (%s, %s)"  % (d_str, base_str, disc_str)

    elif factor_QQ['albert_type'] == 'IV':
        if disc == 1:
            statement = "M_%s (%s)" % (d_str, base_str)
        elif d == 2:
            statement = "Quat (%s, %s)"  % (base_str, disc_str)
        else:
            statement = "Alg_%s (%s, %s)"  % (d_str, base_str, disc_str)

    return statement

def statement_factor_ZZ_maximal(factor_QQ, desc_ZZ, fieldstring):
    # The upcoming value is simply 1 in our current application of this code
    # snippet:
    index = desc_ZZ['index']
    if factor_QQ['albert_type'] == 'I':
        return '%s' % pretty_print_ring(factor_QQ['base_field'], index)
    # TODO: Next line in greater generality
    elif factor_QQ['albert_type'] == 'II' and factor_QQ['disc'] == index and pretty_print_field(factor_QQ['base_field']) == 'QQ':
        return 'M_%s (%s)' % (factor_QQ['dim_sqrt'], pretty_print_ring(factor_QQ['base_field'], index))
    elif factor_QQ['albert_type'] == 'IV' and factor_QQ['disc'] == index and pretty_print_field(factor_QQ['base_field']) == 'QQ':
        return 'M_%s (%s)' % (factor_QQ['dim_sqrt'], pretty_print_ring(factor_QQ['base_field'], index))
    else:
        return 'Max (%s)' % statement_factor_QQ(factor_QQ)

def statement_factors_ZZ_index(factors_QQ, desc_ZZ, fieldstring):
    if len(factors_QQ) == 1:
        statement = "%s" % pretty_print_ring(factors_QQ[0]['base_field'], desc_ZZ['index'])
    else:
        statement = "Sub (End (J_%s) ox QQ, %s)" % (fieldstring, desc_ZZ['index'])
    return statement

def statement_sato_tate_group(desc_dict, genus, fieldstring):
    sato_tate = desc_dict['sato_tate']
    if sato_tate == "":
        sato_tate = "not classified yet"
    return "Sato-Tate group: %s" % sato_tate

def statement_cm(desc_dict, genus, fieldstring):
    factors_QQ = desc_dict['factors_QQ']
    dimsum = 0
    for factor_QQ in factors_QQ:
        if factor_QQ['albert_type'] == 'IV':
            dimsum += (len(factor_QQ['base_field']) - 1) * factor_QQ['dim_sqrt']
    if (dimsum // 2) == genus:
        return "(CM)"
    return ""

def statement_eichler(desc_dict, genus, fieldstring):
    factors_QQ = desc_dict['factors_QQ']
    desc_ZZ = desc_dict['desc_ZZ']
    if len(factors_QQ) == 1 and factors_QQ[0]['albert_type'] != 'I':
        if desc_ZZ['is_eichler'] == 1:
            return "(Eichler)"
    return ""

def statement_gl2(desc_dict, genus, fieldstring):
    factors_QQ = desc_dict['factors_QQ']
    dimsum = 0
    for factor_QQ in factors_QQ:
        dimsum += (len(factor_QQ['base_field']) - 1) * factor_QQ['dim_sqrt']^2
    if dimsum == genus:
        return "of GL_2-type"
    return "not of GL_2-type"

def statement_simple(desc_dict, genus, fieldstring):
    factors_QQ = desc_dict['factors_QQ'] 
    if not len(factors_QQ) == 1:
        factor_QQ = factors_QQ[0]
        if factor_QQ['dim_sqrt'] > 1 and factor_QQ['disc'] != 1:
            return "simple"
    return "not simple"
