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

def pretty_print_over_field(desc):
    desc = sagify_description(desc)
    # Next is the conversion to a dictionary
    # (makes it easier to change as only this part will need modification)
    desc_dict = dict({'factorsQQ' : desc[0], 'factorsRR' : desc[1] , 'factorsZZ' : desc[2]})
    return all_statements(desc_dict, 'F')

def pretty_print_lattice(desc_lat):
    desc_lat = sagify_description(desc_lat)
    # TODO: Apparently such a for-loop is discouraged, but well
    statement = ''
    for desc in desc_lat:
        statement += """Over subfield F with minimal polynomial %s:\n""" % intlist_to_poly(desc[0])
        statement += '\n'
        desc_dict = dict({'factorsQQ' : desc[1], 'factorsRR' : desc[2] , 'factorsZZ' : desc[3]})
        statement += all_statements(desc_dict, 'F')
        statement += '\n'
    return statement

def all_statements(desc_dict, fieldstring):
    statement = ''
    statement += statement_endomorphisms(desc_dict, r'F')
    statement += statement_sato_tate_group(desc_dict)
    statement += statement_gl2(desc_dict)
    statement += statement_simple(desc_dict)
    return statement

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

def intlist_to_poly(s, var = 'x'):
    return str(PolynomialRing(QQ, var)(s))

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
        return 'the number field with defining polynomial %s' % intlist_to_poly(L)

def statement_endomorphisms(desc_dict, fieldstring):
    return "boink\n"
    factorsQQ_number = len(factorsQQ)
    factorsQQ_pretty = [ pretty_print_field(factorQQ[0]) for factorQQ in factorsQQ ]
    statement = """End (J_%s): """ % fieldstring

    # First row: description of endomorphism algebra factors
    statement += """End (J_%s) ox QQ: """ % fieldstring
    # In the case of only one factor we either get a number field or a
    # quaternion algebra:
    if factorsQQ_number == 1:
        # First we deal with the number field case,
        # in which we have set the discriminant to be -1
        if factorsQQ[0][1] == -1:
            # Prettify if labels available, otherwise return defining polynomial:
            statement += """%s""" % factorsQQ_pretty[0]
            # Detect CM by presence of a quartic polynomial:
            if len(factorsQQ[0][0]) == 5 : statement += """ (CM)"""
        # Up next is the case of a matrix ring (trivial disciminant), with
        # labels and full prettification always available:
        elif factorsQQ[0][1] == 1:
            statement += """The 2 x 2 matrix ring over %s""" % factorsQQ_pretty[0]
        # And finally we deal with quaternion algebras over the rationals:
        else:
            statement += """the quaternion algebra over %s of discriminant %s"""\
                % (factorsQQ_pretty[0], factorsQQ[0][1])
    # If there are two factors, then we get two at most quadratic fields:
    else:
        statement += """%s x %s""" % (factorsQQ_pretty[0], factorsQQ_pretty[1])
    # End of second row:
    statement += """\n"""

    # Second row: description of the endomorphism ring as an order in the
    # endomorphism algebra
    # First the case of a maximal order:
    if ring[0] == 1:
        # Single factor:
        if factorsQQ_number == 1:
            # Number field or not:
            if factorsQQ[0][1] == -1:
                # Prettify in quadratic case:
                if len(factorsQQ[0][0]) in [2, 3]:
                    statement += """%s""" % pretty_print_ring(factorsQQ[0][0], 1)
                else:
                    statement += """the maximal order of End (J_%s) ox QQ""" % fieldstring
            else:
                # Use M_2 over integers if this applies:
                if factorsQQ[0][1] == 1 and len(factorsQQ[0][0]) == 2:
                    statement += """M_2 (ZZ)"""
                # TODO: Add flag that indicates whether we are over a PID, in
                # which case we can use the following lines:
                #if factorsQQ[0][2] == 1:
                #    statement += """\(\mathrm{M}_2 (%s)\)"""\
                #        % pretty_print_ring(factorsQQ[0][1], 1)
                else:
                    statement += """a maximal order of End (J_%s) ox QQ""" % fieldstring
        # If there are two factors, then they are both at most quadratic
        # and we can prettify them
        else:
            statement += ' x '.join([ pretty_print_ring(factorQQ[0], 1) for factorQQ in factorsQQ ])
    # Then the case where there is still a single factor:
    elif factorsQQ_number == 1:
        # Number field case:
        if factorsQQ[0][1] == -1:
            # Prettify in quadratic case:
            if len(factorsQQ[0][0]) in [2, 3]:
                statement += """%s""" % pretty_print_ring(factorsQQ[0][0], ring[0])
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
        statement += """an order of index %s in %s""" % (ring[0], ' x '.join([ pretty_print_ring(factorQQ[0], 1) for factorQQ in factorsQQ ]))
    # End of first row:
    statement += """\n"""

    # Third row: description of algebra tensored with RR
    statement += """End (J_%s) ox RR: %s\n""" % (fieldstring, ' x '.join(factorsRR))

    return statement

def statement_sato_tate_group(desc_dict):
    return "boink\n"
    return """Sato-Tate group: %s\n""" % group

def statement_gl2(desc_dict):
    return "boink\n"
    # TODO: Use factorsQQ instead
    # One element in endomorphism algebra, not matrix ring (checked by II with finite discriminant)
    if factorsRR in [ ['RR', 'RR'], ['CC'] ]:
        gl2 = "of GL_2-type"
    else:
        gl2 = "not of GL_2-type"

def statement_simple(desc_dict):
    return "boink\n"
    if len(factorsQQ) == 1 and factorsQQ[0][1] != 1:
        simple = "simple"
    else:
        simple = "not simple"
