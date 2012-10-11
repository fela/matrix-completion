# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 09:19:57 2012

@author: Fela Winkelmolen
"""

import scipy as sp
import scipy.linalg
from matrixgeneration import *
import matrixdecomposition as md

class Tester:
    def __init__(self):
        sp.random.seed(0)
        self._decomposition_args = {}
        self.verbose = True
        self._generate_matrix()
    
    def set_decomposition_options(self, **kwargs):
        self._decomposition_args = kwargs
    
    """ sets the following instance variables:
        * _Y is the input matrix of the decomposition algorithm
        * _Mask is the mask of _Y equal to True in the positions 
          where _Y is observed
        * _TH is the first component of _Y, which we want to recover
        * _GA is the other component"""
    def _generate_matrix(self):
        pass
    
    def comparative_test(self):
        TH, GA = self._TH, self._GA
        Y, Mask, args = self._Y, self._Mask, self._decomposition_args
        A, B = md.matrix_decomposition(Y, Mask, **args)
        if self.verbose == True:
            print GA[1, :]
            print B[1, :]
            print Mask[1, :]
            print TH[1, :]
            print A[1, :]
            
            error = ((TH - A)**2).sum() / TH.size
            print "mean square error:", error
            error2 = ((TH - Y*Mask)**2).sum() / TH.size
            print "mse for naive solution:", error2
            print "improvement:", error2/error, "times"

class SintheticMatrixTester(Tester):
    def __init__(self):
        # default values
        self._size = 30
        self._rank = 3
        self._holes = 0.2 
        self._noise = 5
        Tester.__init__(self)
        
    def _generate_matrix(self):
        size = self._size
        TH = low_rank(size, self._rank) + sp.rand(size, size)/50*self._noise
        GA = spiky(size)
        Mask = selection(size, self._holes)
        Y = (TH + GA) * Mask
        self._Y, self._TH, self._GA, self._Mask = Y, TH, GA, Mask
    
    def no_noise_test(self):
        pass
    
    def default_test(self):
        pass

t = SintheticMatrixTester()
t.comparative_test()