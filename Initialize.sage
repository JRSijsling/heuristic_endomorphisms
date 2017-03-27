"""
 *  Initialization of the Sage part of the package
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
"""

import os
# The following line is a bad solution:
if not '__endodir__' in globals():
    __endodir__ = os.getenv("PWD") + "/"
# The following line is a lazy solution:
magma.chdir(__endodir__)
magma.load('Initialize.m')

# This has to be hidden better to prevent accidental overwriting:
prec = 300
epscomp = 10^(-prec + 30)
epsLLL = 5^(-prec + 7)
epsinv = 2^(-4)
Bound = 48

load(__endodir__ + 'heuristic/Recognition.sage')
load(__endodir__ + 'heuristic/Canonize.sage')
load(__endodir__ + 'heuristic/Decomposition.sage')
load(__endodir__ + 'heuristic/Conversion.sage')
load(__endodir__ + 'Wrapper.sage')
