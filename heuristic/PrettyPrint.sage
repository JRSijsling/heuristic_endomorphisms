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
    return all_statements(desc_dict, genus, fieldstring)

def pretty_print_lattice(desc_lat, genus, fieldstring, varstring):
    desc_lat = sagify_description(desc_lat)
    statements = [ ]
    for desc in desc_lat:
        statement = ''
        field_desc = desc.pop(0)
        statement += """Over subfield F with minimal polynomial %s:\n""" % pretty_print_polynomial_list(field_desc, varstring)
        desc_dict = dict_description(desc)
        statement += all_statements(desc_dict, genus, fieldstring)
        statements.append(statement)
    return '\n\n'.join(statements)

def sagify_description(desc):
    desc_str = str(desc)
    desc_str = desc_str.replace('\n', '').replace('[*', '[').replace('*]', ']').replace(' ','')
    desc_str = desc_str.replace('I', '\'I\'')
    desc_str = desc_str.replace('\'I\'\'I\'', '\'II\'')
    desc_str = desc_str.replace('\'II\'\'I\'', '\'III\'')
    desc_str = desc_str.replace('\'I\'V', '\'IV\'')
    desc_str = desc_str.replace('M_', '\'M_')
    desc_str = desc_str.replace('HH', '\'HH\'')
    desc_str = desc_str.replace('CC', '\'CC\'')
    desc_str = desc_str.replace('RR', '\'RR\'')
    desc_str = desc_str.replace('\'(CC)\'', '(CC)')
    desc_str = desc_str.replace('\'(RR)\'', '(RR)')
    return sage_eval(desc_str)

def dict_factor_QQ(factor_QQ):
    dict_to_fill = dict()
    dict_to_fill['albert_type'] = factor_QQ[0]
    dict_to_fill['base_field'] = factor_QQ[1]
    dict_to_fill['dim_sqrt'] = factor_QQ[2]
    dict_to_fill['discriminant'] = factor_QQ[3]
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
    dict_to_fill['desc_RR'] = dict_desc_RR(desc[1])
    dict_to_fill['desc_ZZ'] = dict_desc_ZZ(desc[2])
    # TODO: Add
    dict_to_fill['sato_tate'] = ""
    return dict_to_fill

def all_statements(desc_dict, genus, fieldstring):
    statements= [ ]
    statements.append(statement_endomorphisms(desc_dict, fieldstring))
    statements.append(statement_sato_tate_group(desc_dict))
    statements.append(statement_gl2(desc_dict, genus) + '; ' + statement_simple(desc_dict))
    return '\n'.join(statements)

def pretty_print_polynomial_list(s, varstring):
    return str(PolynomialRing(QQ, varstring)(s))

def pretty_print_ring(L, f):
    # Only makes a difference for at most quadratic fields
    # FIXME: Generalize, for example to cyclotomic fields
    if len(L) == 2:
        return 'ZZ'
    c,b,a = L
    D = b**2 - 4*a*c
    if f == 1:
        if D % 4 == 0:
            return 'ZZ [sqrt(%s)]'% str(D//4)
        return 'ZZ [(1 + sqrt(%s))/2]' % str(D)
    if D % 4 == 0:
        return 'ZZ [%s sqrt(%s)]' % (str(f), str(D//4))
    if f % 2 == 0:
        return 'ZZ [%s sqrt(%s)]' % (str(f//2), str(D))
    return 'ZZ [(1 + sqrt(%s))/2]' % (str(f), str(D))

def pretty_print_field(L):
    # Only makes a difference for at most quadratic fields
    # FIXME: Generalize, for example to cyclotomic fields
    if len(L) == 2:
        return 'QQ'
    if len(L) == 3:
        c,b,a = L
        D = b**2 - 4*a*c
        return 'QQ (\sqrt(%s))' % D.squarefree_part()
    else:
        return 'the number field with defining polynomial %s' % pretty_print_polynomial_list(L)

def statement_endomorphisms(desc_dict, fieldstring):
    return "boink"
    factors_QQ_number = len(factors_QQ)
    factors_QQ_pretty = [ pretty_print_field(factor_QQ[0]) for factor_QQ in factors_QQ ]
    statement = """End (J_%s): """ % fieldstring

    # First row: description of endomorphism algebra factors
    statement += """End (J_%s) ox QQ: """ % fieldstring
    # In the case of only one factor we either get a number field or a
    # quaternion algebra:
    if factors_QQ_number == 1:
        # First we deal with the number field case,
        # in which we have set the discriminant to be -1
        if factors_QQ[0][1] == -1:
            # Prettify if labels available, otherwise return defining polynomial:
            statement += """%s""" % factors_QQ_pretty[0]
            # Detect CM by presence of a quartic polynomial:
            if len(factors_QQ[0][0]) == 5 : statement += """ (CM)"""
        # Up next is the case of a matrix ring (trivial disciminant), with
        # labels and full prettification always available:
        elif factors_QQ[0][1] == 1:
            statement += """The 2 x 2 matrix ring over %s""" % factors_QQ_pretty[0]
        # And finally we deal with quaternion algebras over the rationals:
        else:
            statement += """the quaternion algebra over %s of discriminant %s"""\
                % (factors_QQ_pretty[0], factors_QQ[0][1])
    # If there are two factors, then we get two at most quadratic fields:
    else:
        statement += """%s x %s""" % (factors_QQ_pretty[0], factors_QQ_pretty[1])
    # End of second row:
    statement += """\n"""

    # Second row: description of the endomorphism ring as an order in the
    # endomorphism algebra
    # First the case of a maximal order:
    if ring[0] == 1:
        # Single factor:
        if factors_QQ_number == 1:
            # Number field or not:
            if factors_QQ[0][1] == -1:
                # Prettify in quadratic case:
                if len(factors_QQ[0][0]) in [2, 3]:
                    statement += """%s""" % pretty_print_ring(factors_QQ[0][0], 1)
                else:
                    statement += """the maximal order of End (J_%s) ox QQ""" % fieldstring
            else:
                # Use M_2 over integers if this applies:
                if factors_QQ[0][1] == 1 and len(factors_QQ[0][0]) == 2:
                    statement += """M_2 (ZZ)"""
                # FIXME: Add flag that indicates whether we are over a PID, in
                # which case we can use the following lines:
                #if factors_QQ[0][2] == 1:
                #    statement += """\(\mathrm{M}_2 (%s)\)"""\
                #        % pretty_print_ring(factors_QQ[0][1], 1)
                else:
                    statement += """a maximal order of End (J_%s) ox QQ""" % fieldstring
        # If there are two factors, then they are both at most quadratic
        # and we can prettify them
        else:
            statement += ' x '.join([ pretty_print_ring(factor_QQ[0], 1) for factor_QQ in factors_QQ ])
    # Then the case where there is still a single factor:
    elif factors_QQ_number == 1:
        # Number field case:
        if factors_QQ[0][1] == -1:
            # Prettify in quadratic case:
            if len(factors_QQ[0][0]) in [2, 3]:
                statement += """%s""" % pretty_print_ring(factors_QQ[0][0], ring[0])
            else:
                statement += """an order of conductor of norm %s in End (J_%s) ox QQ""" % (ring[0], fieldstring)
        # Otherwise mention whether the order is Eichler:
        elif ring[1] == 1:
            statement += """an Eichler order of index %s in End (J_%s) ox QQ""" % (ring[0], fieldstring)
        else:
            statement += """a non-Eichler order of index %s in End (J_{%s}) ox QQ""" % (ring[0], fieldstring)
    # Finally the case of two factors. We can prettify to some extent, since we
    # can describe the maximal order here
    else:
        statement += """an order of index %s in %s""" % (ring[0], ' x '.join([ pretty_print_ring(factor_QQ[0], 1) for factor_QQ in factors_QQ ]))
    # End of first row:
    statement += """\n"""

    # Third row: description of algebra tensored with RR
    statement += """End (J_%s) ox RR: %s\n""" % (fieldstring, ' x '.join(factors_RR))

    return statement

def statement_sato_tate_group(desc_dict):
    sato_tate = desc_dict['sato_tate']
    if sato_tate == "":
        sato_tate = "not classified yet"
    return """Sato-Tate group: %s""" % sato_tate

def statement_cm(desc_dict, genus):
    factor_dicts = desc_dict['factors_QQ']
    dimsum = 0
    for factor_dict in factor_dicts:
        if factor_dict['albert_type'] == 'IV':
            dimsum += (len(factor_dict['base_field']) - 1) * factor_dict['dim_sqrt']
    if dimsum == genus:
        return "(CM)"
    return ""

def statement_gl2(desc_dict, genus):
    factor_dicts = desc_dict['factors_QQ']
    dimsum = 0
    for factor_dict in factor_dicts:
        dimsum += (len(factor_dict['base_field']) - 1) * factor_dict['dim_sqrt']^2
    if dimsum == genus:
        return "of GL_2-type"
    return "not of GL_2-type"

def statement_simple(desc_dict):
    if not len(desc_dict['factors_QQ']) == 1:
        factor_dict = desc_dict['factors_QQ'][0]
        if factor_dict['dim_sqrt'] > 1 and factor_dict['discriminant'] != 1:
            return "simple"
    return "not simple"
