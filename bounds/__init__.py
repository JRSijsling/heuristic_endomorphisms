from sage.all import load

# by using inspect, we can set __endodir__ to the right path, and use heuristic_endomorphisms as a python module
if not '__boundsdir__' in globals():
    import os
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe())[0];
    __boundsdir__ = os.path.dirname(filename) + "/"

from sage.all import *
load(__boundsdir__ + "MagmaInterface.m")
load(__boundsdir__ + "constants.sage");
load(__boundsdir__ + "DiscriminantBound.sage")
load(__boundsdir__ + "TwistPolynomials.sage")
load(__boundsdir__ + "NonQM.sage");
load(__boundsdir__ + "GeometricallyIrreducible.sage");
load(__boundsdir__ + "EndomorphismRankBound.sage")
load(__boundsdir__ + "NonIsogenous.sage")
load(__boundsdir__ + "Genus2Factors.sage")
load(__boundsdir__ + "Test.sage")
