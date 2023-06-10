#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 08:48:14 2022

@author: wentaowu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tmMVFunctions import *

#material properties
H = 0.3          #thickness, m
lamda = 1.5      # thermal conductivity, W/m/K
cp = 750         #specific heat, J/kg/K
rho = 2500       #density, kg/m3

#initial conditions
ACH = 8
Tinit = 24+273.15 # initial temperature of room air and thermal mass

#calculated parameter
alpha = lamda/rho/cp #thermal diffusivity
#ho = 0.13*ACH**0.8    #external convection heat transfer coefficient
ho = 1.9231263 #average
#ho = 0.409958 #min
#ho = 2.857551 #max
Bi = ho*H/lamda      #Biot number

b = transcendental(Bi)
(R,S,W) = tmMVCoeffs(Bi,b)

df_denver = pd.read_csv('WeatherInput_validation.csv')
denverMay1st = df_denver.to_numpy()
Tr = denverMay1st[:,1]
Tmass_max = tmMVTemperature(0.3,H,Tr[0],Bi,alpha,Tr,b,S,W)

#ycoor = np.zeros(31)
#for i in range(1,np.size(ycoor)):
    #ycoor[i] = ycoor[i-1] + H / (np.size(ycoor)-1)

#Tr = np.array([Tinit-273.15,16.1,15.7,15.4,15,14.4,13.9,13.3,13.7,14,14.4,14.6,14.8])
#Tr = Tr + 273.15

#Tinit = 24 + 273.15
#Tr = np.zeros(74)
#Tr[0]  = Tinit
#Tr[1:13] = np.array([16.1,15.7,15.4,15,14.4,13.9,13.3,13.7,14,14.4,14.6,14.8]) + 273.15
#Tr[13:25] = Tinit
#Tr[25:49] = Tr[1:25]
#Tr[49:73] = Tr[1:25]
#Tr[73] = Tr[1]


#TMass = tmMVTemperature(ycoor,H,Tinit,Bi,alpha,Tr,b,S,W)

#hoursPlot = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12])
#plt.plot(hoursPlot, TMass[:,0], label = 'y=0')
#plt.plot(hoursPlot, TMass[:,15], label = 'y=3')

#plt.plot(hoursPlot, TMass[:,30], 'o', label = 'surface')
#plt.ylim(ymin=296.2, ymax = 297.2)
#plt.legend()
#plt.grid()
#plt.xlabel('hours')
#plt.ylabel('Temperature')

#ySurf = np.array([0.3])
#Q = tmMVHeatFlux(ySurf,H,Tinit,Bi,ho,alpha,Tr,b,S,W)
#q = Q[49:]

