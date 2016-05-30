"""
 *  Generates endomorphism data from a colon-separated list.
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

# Usage: Specify indices of defining polynomials and Sato-Tate group here;
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
# The maximum number of retries:
maxrun = 2^5
# The counter for the current run and the line:
run = 0
counter = 0
# Skipping a specific curve if needed:
hell = [ 56306 ]

# Automatic determination of line length of input:
with open(inputfile) as inputstream:
    n = min([ len(line.split(':')) for line in inputstream ])

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
                linestart = ":".join(linesplit[:n])
                # We have to see if there is no new information on the line yet:
                if len(linesplit) == n:
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
                        pol_list = eval(linesplit[fh_index])
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
