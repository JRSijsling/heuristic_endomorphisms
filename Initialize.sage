"""
 *  Initialization
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

if not '__endodir__' in globals():
    raise ImportError("Please set a value for __endodir__")

import os
cur = os.getcwd()

os.chdir(__endodir__)
magma.chdir(__endodir__)

magma.load('~/.magmarc')
magma.AttachSpec('spec')

load('heuristic/Relative.sage')
load('heuristic/PrettyPrint.sage')
load('Wrapper.sage')

import bounds as Bounds

os.chdir(cur)
magma.chdir(cur)
