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
magma.load('~/.magmarc');
magma.AttachSpec('spec');

load(__endodir__ + 'heuristic/Relative.sage')
#load(__endodir__ + 'heuristic/Canonize.sage')
#load(__endodir__ + 'heuristic/Decomposition.sage')
#load(__endodir__ + 'heuristic/Conversion.sage')
load(__endodir__ + 'Wrapper.sage')
