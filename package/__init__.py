# Copyright (C) 2017 Edgar Costa, Jeroen Sijsling
# See LICENSE file for license details.

from sage.all import load

# by using inspect, we can set __endodir__ to the right path, and use heuristic_endomorphisms as a python module
if not '__endodir__' in globals():
    import os
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe())[0];
    __endodir__ = os.path.dirname(filename) + "/"

from sage.all import *
load(__endodir__ + "Initialize.sage");
