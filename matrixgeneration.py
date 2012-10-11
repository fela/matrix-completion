# -*- coding: utf-8 -*-
"""
Created on Thu May 31 11:09:39 2012

@author: Fela Winkelmolen

Some methods to generate test matrixes
"""

import scipy as sp
import scipy.linalg as linalg
import numpy.random as random
#import  load

def delta():
    return sp.finfo(float).eps * 16

def ortonormal(m):
    I = sp.identity(m)
    
    # eigenvectors of a random symmetric m by m matrix
    A = sp.rand(m, m)
    #print '---------------------------'
    #print A.sum()
    S = A * A.T
    _, P = linalg.eigh(S)
    #print P
    #print P.dot(P.T)
    # check if really ortonormal
    if ((P.dot(P.T) - I > delta()).any()):
        return ortonormal(m) # try again
    return P

def low_rank(m, rank):
    diag = sp.concatenate((sp.rand(rank)*5+1, sp.zeros(m-rank)))
    E = sp.diag(diag)
    return ortonormal(m).dot(E).dot(ortonormal(m).T)

def spiky(m, spikes=None, avg=1.0, sigma=None):
    # substitute None with default values
    spikes = spikes or m
    sigma = sigma or avg/6.0
    
    GA = sp.zeros((m, m))
    r = random.randint(0, m, spikes) # row indexes of spikes
    c = random.randint(0, m, spikes) # column indexes of spikes
    GA[[r, c]] = sigma * random.randn(spikes) + avg
    return GA

def selection(m, perc):
    return sp.rand(m, m) >= perc
