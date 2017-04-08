def ReprCurve(X):
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

def ReprEndomorphismData(End):
    return "The endomorphism data of" + ReprCurve(End.X)

def ReprLattice(Lat):
    return "The endomorphism lattice of" + ReprCurve(Lat.X)

def ReprOverField(overfield):
    pre = "The endomorphism structure of" + ReprCurve(overfield.X)
    if overfield.field == "geometric":
        post = " over the algebraic closure of its base field"
    elif overfield.field == "base":
        post = " over its base field"
    else:
        post = " over " + str(overfield.field)
    return pre + post

def ReprDecomposition(decomp):
    return "The decomposition structure of" + ReprCurve(decomp.X)

def PrettyPrintOverField(desc):
    statement = ''
    statement += endo_statement(desc[1], desc[2], desc[3], r'F')
    statement += st_group_statement(desc[4])
    statement += gl2_simple_statement(desc[1], desc[2])
    statement += '\n'
    return "boink"

def PrettyPrintLattice(lat_desc):
    statement = """Smallest field over which all endomorphisms are defined:\nGalois number field K = QQ (a) with defining polynomial %s\n\n""" % intlist_to_poly(lat_desc[0][0])
    for desc in lat_desc:
        statement += """Over subfield F minimal polynomial %s:\n""" % intlist_to_poly(desc[0][0])
        statement += PrettyPrintOverField([ desc[i] for i in [1..len(desc)] ])
    return statement

def intlist_to_poly(s, var = 'x'):
    return str(PolynomialRing(QQ, var)(s))

def ED_sagified(ED):
    return [ [ [ Canonize_Field(factorsQQ[1].sage())[0], factorsQQ[2].sage() ] for factorsQQ in ED[1] ], [ repr(factorRR) for factorRR in ED[2] ], ED[3].sage(), repr(ED[4]) ]

def EDs_sagified(EDs, frep):
    EDs_sage = []
    for ED in EDs:
        ED1 = Canonize_Subfield(ED[1].sage(), frep)
        ED2 = [ [ Canonize_Field(factorsQQ[1].sage())[0], factorsQQ[2].sage() ] for
                factorsQQ in ED[2] ]
        ED3 = [ repr(factorRR) for factorRR in ED[3] ]
        ED4 = ED[4].sage()
        #ED5 = repr(ED[5])
        EDs_sage.append([ED1, ED2, ED3, ED4])
    EDs_sage = sorted(EDs_sage, key = lambda t : len(t[0][0]))
    EDs_sage.reverse()
    return EDs_sage

def ring_pretty(L, f):
    # Only makes a difference for at most quadratic fields
    # FIXME: Generalize, for example to cyclotomic fields
    if len(L) == 2:
        return r'ZZ'
    c,b,a = L
    D = b**2 - 4*a*c
    if f == 1:
        if D % 4 == 0:
            return r'ZZ [sqrt(%s)]'% str(D//4)
        return r'ZZ [(1 + sqrt(%s))/2]' % str(D)
    if D % 4 == 0:
        return r'ZZ [%s sqrt(%s)]' % (str(f), str(D//4))
    if f % 2 == 0:
        return r'ZZ [%s sqrt(%s)]' % (str(f//2), str(D))
    return r'ZZ [(1 + sqrt(%s))/2]' % (str(f), str(D))

def field_pretty(L):
    # Only makes a difference for at most quadratic fields
    # FIXME: Generalize, for example to cyclotomic fields
    if len(L) == 2:
        return r'QQ'
    if len(L) == 3:
        c,b,a = L
        D = b**2 - 4*a*c
        return r'QQ (\sqrt(%s))' % D.squarefree_part()
    else:
        return r'the number field with defining polynomial %s' % intlist_to_poly(L)

def factorsRR_raw_to_pretty(factorsRR):
    return ' x '.join(factorsRR)

def endo_statement(factorsQQ, factorsRR, ring, fieldstring):
    factorsQQ_number = len(factorsQQ)
    factorsQQ_pretty = [ field_pretty(factorQQ[0]) for factorQQ in factorsQQ ]
    statement = """End (J_%s): """ % fieldstring

    # First row: description of the endomorphism ring as an order in the
    # endomorphism algebra
    # First the case of a maximal order:
    if ring[0] == 1:
        # Single factor:
        if factorsQQ_number == 1:
            # Number field or not:
            if factorsQQ[0][1] == -1:
                # Prettify in quadratic case:
                if len(factorsQQ[0][0]) in [2, 3]:
                    statement += """%s""" % ring_pretty(factorsQQ[0][0], 1)
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
                #        % ring_pretty(factorsQQ[0][1], 1)
                else:
                    statement += """a maximal order of End (J_%s) ox QQ""" % fieldstring
        # If there are two factors, then they are both at most quadratic
        # and we can prettify them
        else:
            statement += ' x '.join([ ring_pretty(factorQQ[0], 1) for factorQQ in factorsQQ ])
    # Then the case where there is still a single factor:
    elif factorsQQ_number == 1:
        # Number field case:
        if factorsQQ[0][1] == -1:
            # Prettify in quadratic case:
            if len(factorsQQ[0][0]) in [2, 3]:
                statement += """%s""" % ring_pretty(factorsQQ[0][0], ring[0])
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
        statement += """an order of index %s in %s""" % (ring[0], ' x '.join([ ring_pretty(factorQQ[0], 1) for factorQQ in factorsQQ ]))
    # End of first row:
    statement += """\n"""

    # Second row: description of endomorphism algebra factors
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

    # Third row: description of algebra tensored with RR
    statement += """End (J_%s) ox RR: %s\n""" % (fieldstring, factorsRR_raw_to_pretty(factorsRR))

    return statement

def st_group_statement(group):
    return """Sato-Tate group: %s\n""" % group

def gl2_simple_statement(factorsQQ, factorsRR):
    # TODO: Make this a summation over the RR-type
    if factorsRR in [ ['RR', 'RR'], ['CC'] ]:
        gl2 = "of GL_2-type"
    else:
        gl2 = "not of GL_2-type"
    # One element in endomorphism algebra, not matrix ring:
    if len(factorsQQ) == 1 and factorsQQ[0][1] != 1:
        simple = "simple"
    else:
        simple = "not simple"
    return gl2 + ", " + simple + "\n"
