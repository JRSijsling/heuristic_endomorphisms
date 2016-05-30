"""
 *  Initialization of the Sage part of the package
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

import os
# The following line is a bad solution:
if not __endodir__:
    __endodir__ = os.getenv("PWD")
# The following line is a lazy solution:
magma.chdir(__endodir__)
magma.load('Initialize.m')

# This has to be hidden better to prevent accidental overwriting:
prec = 300
epscomp = 10^(-prec + 30)
epsLLL = 5^(-prec + 7)
epsinv = 2^(-4)
Bound = 48

load(__endodir__ + 'Recognition.sage')
load(__endodir__ + 'Canonize.sage')
load(__endodir__ + 'Decomposition.sage')
load(__endodir__ + 'Conversion.sage')
load(__endodir__ + 'Wrapper.sage')
