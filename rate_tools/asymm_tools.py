"""asymm_tools.py: Functions for calculations needed in strand symmetry modeling."""

#!/usr/bin/python
import threading
import scipy.linalg
from numpy import *
from fp_check import *
from RptMatrixMod import *
#from semiphore import *

############################################################
# MatCalcExcep: An exception to throw when a calculation isn't working right.
class MatCalcExcep(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return str(self.msg)

############################################################
# Usful matrix functions
def matrix2list(M):
    """Convert a matrix two a row-order list"""
    return list(M.flat)

def list2matrix(L, num_rows = -1):
    """Convery a list representing a row-major ordering of a matrix to a matrix.
       Assumed to be square if not otherwise specified."""
    M = matrix(L)
    return M.reshape(sqrt(len(L)) if num_rows == -1 else num_rows, -1)

def matrixSum(L, dtype = None):
    """Return the sum of a list of matrices"""
    S = zeros(L[0].shape, dtype if dtype else L[0].dtype)
    for M in L:
        S += M
    return S

def matrixAverage(L, dtype = None):
    """Take a list of same-sized matrices and return the average matrix"""
    return matrixSum(L, dtype) / len(L)

def isP(M):
    """Return True if M is a legitimate discrete transition matrix.
    (That is: every element is a probability, and row-sums are 1.)"""
    return isOne(M.sum(1)).all() and gteZero(M).all() and gteZero(-1*M + 1).all()

def isQ(M):
    """Return True if M is a legitimate continuous MC matrix matrix.
    (That is, every off-diagonal element is non-negative, and row-sums are 0.)
    """
    return isZero(M.sum(1)).all() and gteZero(-1*diag(M)).all() and gteZero(M - diag(M)*identity(4)).all()


############################################################
# Functions for computing the different realted matrices
def compute_P(C, asymm = True):
    """Given a C matrix compute the P matrix.
    * C: The count-matrix.
    * asymm: boolean reflecting whether the model should be strand asymmetric.

    Raises MatCalcExcept if:
    * Matrix has a zero diagonal element.
    * Other computation error.
    """
    if not C.diagonal().all():
        raise MatCalcExcep("Count matrix has a zero diagonal (raised by compute_P).")

    if asymm:
        P = vstack([C[i,:]/float64(sum(C[i,:])) if sum(C[i,:]) > 0 else zeros((1,4)) for i in range(4)])
    else:
        P = zeros((4,4), float64)

        lambda1 = float64(C[0,0] + C[3,3] + C[0,1] + C[3,2] + C[0,2] + C[3,1] + C[0,3] + C[3,0])
        lambda2 = float64(C[1,0] + C[2,3] + C[1,1] + C[2,2] + C[1,2] + C[2,1] + C[1,3] + C[2,0])

        P[0,0] = (C[0,0] + C[3,3]) / lambda1
        P[0,1] = (C[0,1] + C[3,2]) / lambda1
        P[0,2] = (C[0,2] + C[3,1]) / lambda1
        P[0,3] = (C[0,3] + C[3,0]) / lambda1

        P[1,0] = (C[1,0] + C[2,3]) / lambda2
        P[1,1] = (C[1,1] + C[2,2]) / lambda2
        P[1,2] = (C[1,2] + C[2,1]) / lambda2
        P[1,3] = (C[1,3] + C[2,0]) / lambda2

        P[2,0] = P[1,3]
        P[2,1] = P[1,2]
        P[2,2] = P[1,1]
        P[2,3] = P[1,0]

        P[3,0] = P[0,3]
        P[3,1] = P[0,2]
        P[3,2] = P[0,1]
        P[3,3] = P[0,0]
        
    return P    

def compute_Rt(P):
    """Given a P matrix, create the matrix Rt = logm(P) -- the continuous 
    analog to the P discrete transition matrix."""
    try:
        Rt = matrix(scipy.linalg.logm(P, disp = False)[0])
    except Exception as E:
        raise MatCalcExcep(str(E))

    return makeReal(Rt)

def compute_d(M):
    """Compute the average rate of change of matrix M, where M is the matrix
    for either a discrete or continuous Markov Chain."""
    if isP(M):
        return -0.25*(linalg.slogdet(M)[1])
    elif isQ(M):
        return -0.25*mean(diagonal(M))
    else:
        raise MatCalcExcep("compute_d given a matrix that was not a correct probability or rate matrix")

def compute_q(M, d = None):
    """Compute the q matrix from the P or Rt matrix.  If d not provide, calculates it."""
    if d == None:
        d = compute_d(M)

    if isP(M):
        Rt = compute_Rt(M)
    elif isQ(M):
        Rt = M
    else:
        raise MatCalcExcep("compute_q given a matrix that was not a correct probability or rate matrix")

    if isZero(d) or d == inf:
        raise MatCalcExcep("Computing q for a zero- or inf-distance matrix (compute_q)")

    return Rt / d
                               
