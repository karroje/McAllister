"""
fp_check.py:
Functions for comparing floating-point values, using a tolerance (default value 0.0001) to
define euqality.  (e.g. x==y if abs(x - y) < tolerance0.

For each function there is a vectorized version, which applies the function to each element 
of a matric, returning a matric of the results.  (Having learned more of NumPy since writing this,
its not clear to me these are necessary -- seems to be automatic.)
"""

import sys
import math
from numpy import matrix
from numpy import vectorize

_tolerance = 0.00001
def set_tolerance(t):
    """Change the tolerance constant (default = 0.00001)"""
    global _tolerance
    _tolerance = t


def isZero(a):
    """Determine if a floating point is 0 within a floating-point _tolerance factor"""
    return abs(a) < _tolerance

isZeroV = vectorize(isZero)



def isOne(a):
    """Check if a floating point is 1 within a floating-point _tolerance factor"""
    return abs(a-1) < _tolerance

isOneV = vectorize(isOne)



def gtZero(a):
    """Determine if a is greater than 0 with a floating-point _tolerance factor"""
    return a >= _tolerance

gtZeroV = vectorize(gtZero)



def gteZero(a):
    """Determine if a is >= 0 with a a floating-point _tolerance factor"""
    return a >= -_tolerance

gteZeroV = vectorize(gteZero)



def isReal(a):
    """Determine if a floating point is a real number within a floaing-point _tolerance factor"""
    return abs(a.imag) < _tolerance

isRealV = vectorize(isReal)


def makeReal(a):
    """Convert a compex number to a real by truncating the imaginary portion."""
    return a.real

makeRealV = vectorize(makeReal)


def isEqual(a, b):
    """Determine if two floating points are equal within a floating-point _tolerance factor"""
    return abs(a-b) < _tolerance

def isEqualV(a,b):
    """Vectorized version of isEqual"""
    return _isZero(a-b)

