# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 09:19:57 2012

@author: Fela Winkelmolen

Tests using the data from http://archive.ics.uci.edu/ml/datasets/SPECTF+Heart
"""

import scipy as sp
import scipy.linalg
from matrixgeneration import *
import matrixdecomposition as md


#############################################################
class Experiment:
    def __init__(self, **opt):
        # default values
        self.filename = opt.pop('filename', None) or 'SMALL'
        self.seed = opt.pop('seed', None) or 0
        self.holes = opt.pop('holes', None)
        # do not overwrite if holes is set to .0!
        if self.holes == None:
            self.holes = 0.2 # holes percentage
        self.completion = opt.pop('completion', None) or 'matrix'
        
        self.options = opt
        self.generate_matrix()

    def run(self):
        self.generate_mask()
        TH, GA = self.TH, self.GA
        Y, Mask = self.Y, self.mask
        A = self.matrix_completion()
        
        n = (1-Mask).sum()
        diff = (TH-A)*(1-Mask)
        sqerr = (diff**2).sum()
        if n != 0:
            sqerr /= n
        return sqrt(sqerr)
    
    def matrix_completion(self):
        if self.completion == 'matrix':
            res, _ = md.matrix_decomposition(self.Y, self.mask, **self.options)
        elif self.completion == 'mean':
            means = self.Y.mean(0) # means along the column axis
            # calculate the true means along the column axis
            # using only the values in the mask
            rows, cols = self.Y.shape
            means = sp.zeros(cols)
            # calculate the true means along the column axis
            # using only the values in the mask
            for c in range(cols):
                sum = 0.0
                n = 0
                for r in range(rows):
                    if self.mask[r][c]:
                        sum += self.Y[r][c]
                        n += 1
                if n == 0:
                    means[c] = 0
                else:
                    means[c] = sum/n
            
            res = double(self.Y)
            for r in range(rows):
                for c in range(cols):
                    if not self.mask[r][c]:
                        res[r][c] = means[c]
        return res
        
    
    # called once after initialization
    def generate_matrix(self):
        filename = 'data/' + self.filename + '.train'
        self.TH, _ = self.load_data(filename)
        self.GA = self.TH * 0
        self.Y = self.TH + self.GA
    
    def generate_mask(self, holes=None):
        if holes != None:
            self.holes = holes
        rows, cols = self.GA.shape
        if self.seed != None:
            seed(self.seed)
        self.mask = rand(rows,cols) > self.holes
    
    @staticmethod
    def load_data(filename):
        M = [[int(i) for i in line.split(',')] for line in open(filename)]
        M = sp.array(M)
        y = M[:, 0]
        X = M[:, range(1, M.shape[1])]
        return (X, y)

def holes_experiment(**opt):
    n = opt.pop('steps', None) or 5
    runs = opt.pop('runs', None) or 1
    label = opt.pop('label', None)
    #opt['mu_d'] = 0
    y = sp.array(range(n+1)) / float(n)
    x = []
    for holes in y:
        print '.',
        acc = []
        e = Experiment(**opt)
        for i in range(runs):
            e.generate_mask(holes)
            acc.append(e.run())
        x.append(sp.array(acc).mean())
    plot(y, x, label=label)
    legend(loc=0)

def param_experiment(param_name, params, **opt):
    label = opt.pop('label', None) or param_name
    scale = opt.pop('scale', None) or 'linear'
    x = []
    for p in params:
        print '.',
        opt[param_name] = p
        e = Experiment(**opt)
        x.append(e.run())
    xscale(scale)
    plot(params, x, label=label)
    legend(loc=0)



# exponential range
def exp_range(minval=0.001, maxval=100, steps=10):
    min_exp = sp.log(minval)
    max_exp = sp.log(maxval)
    return sp.exp(sp.linspace(min_exp, max_exp, num=steps))

# test for different completion percentage
# more runs are made to get an estimate of the variance
# and completion using the mean of the column is used for comparison
def experiment1():
    for s in range(5):
        holes_experiment(steps=10, alpha=100000, completion='matrix', seed=s)
    
    #holes_experiment(steps=10, alpha=100000, mu_d=1, completion='matrix', label='mu_d=1')
        
    holes_experiment(steps=20, runs=5, alpha=100000, completion='mean', seed=0, label='mean')

# test different values of mu_d
def experiment2():
    params = exp_range(0.00001, 100, 30)
    param_experiment('mu_d', params, alpha=100000, label='0.2')
    figure()
    params = exp_range(0.00001, 100, 30)
    param_experiment('mu_d', params, alpha=100000, holes=0.6, label='0.6')

# test different values of lambda_d
def experiment3():
    params = exp_range(0.005, 0.2,)
    param_experiment('lambda_d', params, alpha=100000, label='0.2')
    figure()
    param_experiment('lambda_d', params, alpha=100000, holes=0.6, label='0.2')
    