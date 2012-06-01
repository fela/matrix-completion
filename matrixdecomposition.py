# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 12:17:57 2012

@author: Fela Winkelmolen
"""


import scipy as sp
import scipy.linalg

def matrix_decomposition(Y, lambda_d, mu_d, alpha, X=None):
    t = 1.0
    # dimensions are correct supposing X being the Identity matrix
    TH = GA = TH_last = GA_last = Y * 0.0 # theta and gamma
    for k in range(0, 3000):
        t_last = t
        #print (t_last - 1)/t
        t = ( 1.0 + sp.sqrt(1 + 4 * t**2) )/ 2.0
        #print (t_last - 1)/t
        Z = TH + (t_last - 1)/t*(TH-TH_last)
        N = GA + (t_last - 1)/t*(GA-GA_last)
        lambd = 0.5 # TODO: check
        f = ((Z+N) - Y) * 0.5
        first = Z - f
        second = N - f
        TH_last, GA_last = TH, GA
        TH, GA = prox_g([first, second], lambd, lambda_d, mu_d, alpha)
        if (abs(TH - TH_last) + abs(GA - GA_last)).max() < 1e-12:
            print k, "iterations"
            break
    return [TH, GA]

def prox_g(grad, lambd, lambda_d, mu_d, alpha):
    N, Z = grad
    
    X = N
    P = Q = N*0.0
    for i in range(0, 3): # TODO: how many iterations?
        Y = prox1(X+P, lambda_d, lambd)
        P = X + P - Y
        X_last = X
        X = prox2(Y+Q, alpha)
        Q = Y + Q - X
    V = X
    
    # soft thresholding
    W = soft_threshold(Z, mu_d*lambd)
    return [V, W]

# projection of nuclear norm
def prox1(X, lambda_d, lambd):
    U, s, Vh = sp.linalg.svd(X)
    d1, d2 = X.shape
    E = sp.linalg.diagsvd(s, d1, d2) # sigma 
    S = soft_threshold(E, lambda_d*lambd)
    return U.dot(S.dot(Vh))

# projection in Q
def prox2(X, alpha):
    limit = alpha / sp.sqrt(sp.array(X.shape).prod())
    #X = sp.minimum(X, limit)
    #X = sp.maximum(X, -limit)
    return X

def soft_threshold(X, s):
    return (X-s)*(X>s) + (X+s)*(X<-s)


#Y = sp.array([[1, 2, 3], [2, 100, 6]])
#A, B = matrix_decomposition(Y, 0.1, 0.08, 1000)
#print A
#print B
#print A+B