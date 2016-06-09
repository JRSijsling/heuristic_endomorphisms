Description
-----------

This repository contains a mix of Magma, Pari and Sage code for calculating a heuristic (and usually correct) approximation of the endomorphism algebras and rings of Jacobian varieties of hyperelliptic curves. For now the code assumes the curves involved to be defined over QQ, and most of it additionally assumes that we are in genus 2.

Installation
------------

An installation of both Magma and Sage is required to run this code. It can be loaded inside Sage by going to the package/ directory and typing

load('Initialize.sage')

Alternatively, add the following lines to the Sage initialization file (typically found in ~/.sage/init.sage) to enable the functionality on startup anywhere in the system:

\_\_endodir\_\_ = '~/[PATH]/package/'  
load(\_\_endodir\_\_ + 'Initialize.sage')

where [PATH] is the path to the cloned or copied repository.

A bug fix
---------

It is highly recommended (though not absolutely necessary) to fix a Magma bug before using this package. In old version the file magma/package/Algebra/AlgQuat/interface.m had the following as line 145:

c := [Trace(theta), Norm(theta)];

This should be replaced by

cpol := MinimalPolynomial(theta);  
assert Degree(cpol) eq 2;  
c := [Coefficient(cpol,1), Coefficient(cpol, 0)];

Usage
-----

A sample run is given in package/Examples.sage.
