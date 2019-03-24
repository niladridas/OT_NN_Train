#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 15:51:44 2019

@author: vish0908
"""

import numpy as np
import scipy as sp
# q = [alpha,beta,kappa]
# xM: Mean of state x
# xC: Covariance of statex
# vM: Mean of process noise
# vC: Covariance of process noise
# nM: Mean of measurement noise
# nC: Covariance of measurement noise

def sigma(xM,xC,vC,nC,q,t):
    Pa1 = np.concatenate((xC,np.zeros(np.shape(xC)[1],np.shape(vC)[2]),np.zeros(np.shape(xC)[1],np.shape(nC)[2])), axis = 1)
    Pa2 = np.concatenate((np.zeros(np.shape(vC)[1],np.shape(xC)[2]),vC,np.zeros(np.shape(vC)[1],np.shape(nC)[2])), axis = 1)
    Pa3 = np.concatenate((np.zeros(np.shape(nC)[1],np.shape(xC)[2]),np.zeros(np.shape(nC)[1],np.shape(vC)[2]),nC), axis = 1)
    Pa  = np.concatenate((Pa1,Pa2,Pa3)) # Pa is the diagonal-like matrix from the algorithm
    
    x0a_hat = np.concatenate((xM,vM,nM)) # Augmented state
    
    L = np.shape(xM)[1]
    lam = (q[1]**2)*(L + q[3]) - L
    chiA = np.concatenate((x0a_hat,x0a_hat+sp.linalg.sqrtm((L + lam)*Pa), x0a_hat - sp.linalg.sqrtm((L+lam)*Pa)), axis = 1)
        
    if t = 1:
        wm = lam/(lam + L)
        wc = lam/(lam + L) + (1-q[1]**2 + q[2])
    else:
        wm = 1/(2*(lam + L))
        wc = wm
    sigW = {"sigma": chiA, "wm": wm, "wc": wc}
    return sigW