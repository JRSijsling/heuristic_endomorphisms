"""
 *  Initialization of the Sage part of the package
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

if not '__endodir__' in globals():
    import os
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe())[0];
    __endodir__ = os.path.dirname(filename) + "/"

magma.load('~/.magmarc')
magma.AttachSpec(__endodir__ + 'spec')

load(__endodir__ + 'heuristic/Relative.sage')
#load(__endodir__ + 'heuristic/Conversion.sage')
#load(__endodir__ + 'heuristic/Decomposition.sage')
load(__endodir__ + 'Wrapper.sage')

# This is very suboptimal; need relative import
load(__endodir__ + 'bounds/__init__.py')
#import bounds as Bounds
