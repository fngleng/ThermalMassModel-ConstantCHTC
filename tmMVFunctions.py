#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 18:36:47 2022

@author: wentaowu
"""

import numpy as np
import matplotlib.pyplot as plt
import pychebfun as cf

def transcendental(Bi):
    def transEq(x):
        return (-x*np.sin(x)+Bi*np.cos(x))
    f_cheb = cf.Chebfun.from_function(transEq, domain = (0,500))
    rts = f_cheb.roots()
    rts.sort()
    xcorr = np.linspace(0,10,100)
    plt.plot(xcorr, transEq(xcorr), label = 'function')
    plt.plot(rts[0:3], transEq(rts[0:3]), 'o', label = 'roots')
    plt.ylim(ymin=-10, ymax = 10)
    plt.legend()
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    return rts
    
def tmMVCoeffs(Bi,b):
    N = np.size(b)
    R = np.zeros(N)
    S = np.zeros(N)
    W = np.zeros(N)
    for i in range(0,N):
        R[i] = 0.5 * (Bi**2 + b[i]**2 + Bi) / (Bi**2 + b[i]**2)
        S[i] = Bi * np.sin(b[i]) / (b[i] * R[i])
        W[i] = Bi * (np.cos(b[i]) - np.sin(b[i])/b[i]) / (R[i] * b[i]**2)
    return (R, S, W)

def tmMVTemperature(ycoor,H,Tinit,Bi,alpha,Tr,b,S,W):
    
    K = np.size(ycoor) #number of locations to calculate T
    hours = np.size(Tr)    #hours in total
    thetaR= (Tr - Tinit) / Tinit
    
    N = np.size(b)     #Fourier items
    eta = ycoor / H

    theta = np.zeros([hours,K])
    TMass = np.zeros([hours,K])
    TMass[0,:] = Tinit
    
    for l in range(0,K):
        
        for M in range(1,hours):
            temp = 0.0
            tauM = M * 3600 * alpha / (H**2)
            tau1 = 1 * 3600 * alpha / (H**2)
            
            for n in range(0,N):
                exp1 = np.exp(-(b[n]**2) * (tauM-tau1))
                expM = np.exp(-(b[n]**2) * tauM)
                temp1 = thetaR[1] * (exp1 - expM)
                temp2 = thetaR[1] * expM
        
                for m in range(1,M):
                    taum = m * 3600 * alpha / (H**2)
                    taum1 = (m+1) * 3600 * alpha / (H**2)
                    expm0 = np.exp(-(b[n]**2) * (tauM-taum))
                    expm1 = np.exp(-(b[n]**2) * (tauM-taum1))
                   
                    temp1 = temp1 + thetaR[m+1] * (expm1 - expm0)
                    temp2 = temp2 + (thetaR[m+1]-thetaR[m]) * expm0
                
                temp = temp + np.cos(b[n]*eta) * (S[n] / b[n]**2 * temp1 
                                                     - W[n] * temp2)
            theta[M,l] = temp + Bi * thetaR[M] * 0.5 * (eta**2 - 1)
            TMass[M,l] = Tinit + Tinit *  theta[M,l]
    return TMass

def tmMVHeatFlux(ySurf,H,Tinit,Bi,ho,alpha,Tr,b,S,W):
    
    hours = np.size(Tr)
    Q = np.zeros(hours)
    TMass = tmMVTemperature(ySurf,H,Tinit,Bi,alpha,Tr,b,S,W)
    Q = ho * (Tr - TMass[:,0])
    return Q

    
    
    
    
    
    