"""
 *  Generates endomorphism data from a colon-separated list.
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

# Adds endomorphism data to a file of colon-separated lines

# Specify index of defining polynomials:
fh_index = 3
# Precision (not below 200 please):
prec = 100

import os, shutil

# Specify input and output:
inputfile = 'gce_1000000-input.txt'
outputfile = 'gce_1000000-output.txt'

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
                End = EndomorphismData(X, prec = prec, have_oldenburg = True)
                P = End.period_matrix()
                #Lat_desc = End.lattice().descriptions()
                #outputstream.write(linestart + ':' + str(Lat_desc[1][len(Lat_desc[1])]) + '\n')
            except:
                print "Error"
                outputstream.write(linestart + ':' + 'Error' + '\n')
