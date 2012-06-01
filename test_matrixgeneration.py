# -*- coding: utf-8 -*-
"""
Created on Thu May 31 11:18:10 2012

@author: Fela Winkelmolen
"""

import unittest
import matrix_generation as mg

import scipy as sp
import numpy.random as random
import numpy.linalg

class TestMatrixGeneration(unittest.TestCase):
    def test_ortonormal(self):
        n = 15
        I = sp.identity(n)
        for _ in range(0, 100):
            M = mg.ortonormal(n)
            self.assertTrue( (M.dot(M.T) - I <= mg.delta()).all() )
    
    def test_low_rank(self):
        for _ in range(0, 100):
            rank = random.randint(3, 8)
            M = mg.low_rank(15, rank)
            actual_rank = numpy.linalg.matrix_rank(M, mg.delta() * 4)
            self.assertEqual(actual_rank, rank)

if __name__ == '__main__':
    unittest.main()
