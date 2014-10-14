"""rpt_matrix.py: A RptMatrix class for storing needed information about individual repeats along with 
their C matrix.  Will store these in psm files (pickled substitution matrix)"""

from RMfileReader import *
import os.path
import argparse
import pickle
import glob
import re
from numpy import zeros

base2index = {'A':0,'C':1,'G':2,'T':3}                 # Dictionary mapping bases to their index position.
baseComplement = {x:y for x,y in zip("ACGT","TGCA")}   # Dictionary mapping bases to their complements. 
baseSet = {'A','C','G','T'}                            # Set of legal bases
def reverse_complement(s):
    return "".join([baseComplement[c] if c in baseComplement else c for c in s[::-1]])

class RptMatrix:
    def __init__(self, r, M):
        """
        Class for holding both a repeat object and its calculated substitutinon matrix.  Has the following atrributes:
        * ancestor_coords: A pair specifying its [start,finish) coordinates relative to the ancestral sequence.
        * start, finish: the [start,finish) interval on the chromosome, relative to the positive strand.
        * chr: The chromosome name (e.g. "chr1").
        * chr: The containing strand (+ or C).
        * rep_name: The name of the general familiy of repeat types.
        * class_name: The name of the specific type of repeat.
        * M: The 4x4 matrix count of aligned bases, representing possibilities in alphabetical order.
        """
        if r is None:   # Allows us to create an "empty" object"
            assert M is None   # Empty objects shouldn't be given a matrix
            self.ancestor_coords = (None, None)
            self.start, self.finish = None, None
            self.chr = None
            self.strand = None
            self.rep_name = None
            self.class_name = None
            self.M = zeros((4,4))
        else:
            self.ancestor_coords = r.ancestor_coords()
            self.start, self.finish = r.modern_coords()    # Coordinates
            self.chr = r.modern_chr()
            self.strand = r.modern_strand()
            self.rep_name = r.rep_name()
            self.class_name = r.class_name()        
            self.M = M

    def __getitem__(self, index):
        return self.M[index[0], index[1]]

def load_psm(file):
    """Load in a list of RptMatrix objects from a psm file."""
    return pickle.load(open(file, "rb"))

def rptsInPartition(rpt_list, start, finish):
    """Return a sublist of all repeats contained within the [start, finish) interval"""
    i = 0
    while i < len(rpt_list) and rpt_list[i].start < start:
        i += 1
    j = i+1
    while j < len(rpt_list) and rpt_list[j].finish <= finish:
        j = j+1
    return rpt_list[i:j]    






