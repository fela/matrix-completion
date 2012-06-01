# -*- coding: utf-8 -*-
"""
Created on Thu May 31 11:09:39 2012

@author: fela
"""

import scipy as sp
import scipy.linalg as linalg

def delta():
    return sp.finfo(float).eps * 4

def ortonormal(m):
    I = sp.identity(m)
    
    # eigenvectors of a random symmetric m by m matrix
    A = sp.rand(m, m)
    S = A * A.T
    _, P = linalg.eigh(S)
    # check if really ortonormal
    if ((P.dot(P.T) - I > delta()).any()):
        return ortonormal(m) # try again
    return P

def low_rank(m, rank):
    diag = sp.concatenate((sp.rand(rank)*5+1, sp.zeros(m-rank)))
    E = sp.diag(diag)
    return ortonormal(m).dot(E).dot(ortonormal(m).T)
