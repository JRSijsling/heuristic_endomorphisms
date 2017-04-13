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

# Defining polynomials have to be provided in pairs, defined by strings in x or
# by lists of integers. These polynomials (and the conjectural Sato-Tate group,
# if provided) need to be at a consistent index in the provided lines.

# Specify indices of defining polynomials and Sato-Tate group here;
# making the latter negative ignores the final check.
fh_index = 1
st_index = -1
# Precision (not below 200 please):
prec = 300

import os, shutil

# Specify input and output:
inputfile = 'database_input.txt'
intermediatefile = 'database_intermediate.txt'
outputfile = 'database_output.txt'

# Ambient ring needed for substitution:
R.<x> = PolynomialRing(QQ)
# Safeguards against Van Wamelen's infinite loops and other errors:
# First substitutions (do not touch or an infinite loop will result):
subst_list = [ x + 1, x ]
# Bound on height of random substitutions after these have been tried:
B = 3
# The maximum number of retries for difficult curves:
maxrun = 2^5
# Indices of particularly nasty curves (typically those for which the period
# calculation leads into an infinite loop):
hell = [ 56306 ]
# The counter for the current run and the line:
run = 0
counter = 0
done_list = [ ]

stop = False
exhaust = False
while not stop:
    run += 1
    counter = 0
    exhaust = len(subst_list) == 0
    if not exhaust:
        subst = subst_list.pop()
    else:
        while True:
            num = R([QQ.random_element(B) for i in range (2)])
            den = R([QQ.random_element(B) for i in range (2)])
            if num != 0 and den != 0:
                break
            # We do not have to guard against the constant case, since that
            # bugs in Magma already and will get caught below.
        subst = num/den
    stop = True
    print "Run:", run
    print "Substitution:", subst
    with open(inputfile) as inputstream:
        with open(outputfile, 'w') as outputstream:
            for line in inputstream:
                counter += 1
                linestrip = line.rstrip()
                linesplit = linestrip.split(':')
                linestart = linestrip
                # We have to see if there is no new information on the line yet:
                if not counter in done_list:
                    # In the USp(4) case we know everything:
                    if linesplit[st_index] == "USp(4)":
                        outputstream.write(linestart + ':' +
                                # Modify this if the data display changes
                                #"[[[[0,1],[0]],[[[0,1],-1]],['RR'],[1,-1],'USp(4)']]:[[0,1],[0]]:[]"
                                "[[[[0,1],[0]],[[[0,1],-1]],['RR'],[1,-1],'USp(4)']]:[]:[]"
                                + '\n')
                    # Avoiding a nasty infinite loop:
                    elif (subst == x) and (counter in hell):
                        outputstream.write(line)
                    else:
                        print counter
                        pol_list = eval(linesplit[fh_index].replace('^', '**'))
                        f = R(pol_list[0])
                        h = R(pol_list[1])
                        den = subst.denominator()
                        f = R(den^6 * f(subst))
                        h = R(den^3 * h(subst))
                        try:
                            End = EndomorphismData(f, h, prec = prec)
                            Lat = End.lattice()
                            Lat_str = End._lat_
                            Dec = End.decomposition()
                            ECs = Dec.factors()
                            ECs_str = Dec._ECs_rep_
                            if ECs_str:
                                SplFoD_str = Dec._spl_fod_
                            else:
                                # TODO: Check compatibility of next line with LMFDB
                                SplFoD_str = []
                            # Check Sato-Tate and write if a match occurs:
                            if st_index < 0 or Lat_str[-1][4] == linesplit[st_index]:
                                outputstream.write(linestart
                                        + ':' + repr(Lat_str).replace('\n', '').replace(' ', '')
                                        + ':' + repr(SplFoD_str).replace('\n', '').replace(' ', '')
                                        + ':' + repr(ECs_str).replace('\n', '').replace(' ', '')
                                        + '\n')
                            else:
                                outputstream.write(line)
                            done_list.append(counter)
                        except:
                            # In case of an error postpone until next time:
                            outputstream.write(line)
                            stop = False
                else:
                    # Skip the line if it has been calculated already:
                    outputstream.write(line)
    inputfile = intermediatefile
    shutil.copyfile(outputfile, inputfile)
    if run >= maxrun:
        stop = True
os.remove(intermediatefile)
