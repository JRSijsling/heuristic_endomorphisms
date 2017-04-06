# Description

This repository contains a mix of Magma, Pari and Sage code for calculating a heuristic (and usually correct) approximation of the endomorphism algebras and rings of Jacobian varieties of hyperelliptic curves.
For now, the code assumes the curves involved to be defined over QQ, and most of it additionally assumes that we are in genus 2.

# Installation

An installation of both Magma and Sage is required to run this code. 

To install the package pick a directory where you would like to place a copy of the package and do:
```
git clone https://github.com/JRSijsling/heuristic_endomorphisms.git
```

# Loading the package
The `heuristic_endomorphisms` package can be loaded in Sage in three ways.
Let `[PATH]` denote the path for the directory into which you have cloned or copied the repository.

### 1. As a python module
```
sys.path.append('[PATH]')
from heuristic_endomorphisms import *
``` 
to make it more convenient you could also add the line
```
sys.path.append('[PATH]')
```
the Sage initialization file (typically found in  `~/.sage/init.sage`).

### 2. Using `load()`
By going to `[PATH]/heuristic_endomorphisms` directory and typing
```
load('Initialize.sage')
```

### 3. Loaded at startup 
If you prefer, you can have Sage to load every time you start Sage by adding the following lines to the Sage initialization file (typically found in ~/.sage/init.sage):
```
__endodir__ = '[PATH]/heuristic_endomorphisms'  
load(__endodir__ + 'Initialize.sage')
```

alternatively

```
sys.path.append('[PATH]')
from heuristic_endomorphisms import *
```

Note that this will also startup Magma automatically.

# Usage 

A sample run is given in `Examples.sage`.

# A bug fix


It is highly recommended to fix a Magma bug before using this package. In old version the file `magma/package/Algebra/AlgQuat/interface.m` had the following as line 145:
```
c := [Trace(theta), Norm(theta)];
```
This should be replaced by
```
cpol := MinimalPolynomial(theta);  
assert Degree(cpol) eq 2;  
c := [Coefficient(cpol,1), Coefficient(cpol, 0)];
```
