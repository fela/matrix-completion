# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:18:35 2012

@author: Fela Winkelmolen
"""

import unittest
import matrixdecomposition as md

import scipy as sp

class TestMatrixGeneration(unittest.TestCase):
    def test_decomposition_id(self):
        # testing using the identity map
        TH = sp.array([[1, 2, 3], [2, 4, 6]])
        GA = sp.array([[0, 0, 0], [0, 100, 0]])
        
        Y = TH + GA
        
        A, B = md.matrix_decomposition(Y, lambda_d=0.1, mu_d=0.08, alpha=1000)
        
        self.assertLess(abs(TH-A).sum(), 1.0, "")
        self.assertLess(abs(TH-A).max(), 0.2, "")
        
        self.assertLess(abs(GA-A).sum() < 1.0, "")


if __name__ == '__main__':
    unittest.main()
