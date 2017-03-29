"""
 *  Generates endomorphism data from a colon-separated list.
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
"""

# Adds endomorphism data to a file of colon-separated lines

# Specify index of defining polynomials:
fh_index = 2
# Precision (not below 200 please):
prec = 300

import os, shutil

# Specify input and output:
inputfile = 'gce_genus3_hyperelliptic_cond-input.txt'
outputfile = 'gce_genus3_hyperelliptic_cond-output.txt'

# Ambient ring needed for substitution:
R.<x> = PolynomialRing(QQ)
# The counter for the line:
counter = 0
done_list = [ ]

with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        for line in inputstream:
            counter += 1
            linestrip = line.rstrip()
            linesplit = linestrip.split(':')
            linestart = linestrip
            print counter
            pol_list = eval(linesplit[fh_index].replace('^', '**'))
            f = R(pol_list[0])
            h = R(pol_list[1])
            X = HyperellipticCurve(f, h)
            try:
                End = EndomorphismData(X, prec = prec)
                Lat = End.lattice()
                Lat_str = repr(End._lat_[1])
                outputstream.write(linestart + ':' + Lat_str.replace('\n', '').replace(' ', '') + '\n')
            except:
                print "Error"
                outputstream.write(linestart + ':' + 'Error' + '\n')

