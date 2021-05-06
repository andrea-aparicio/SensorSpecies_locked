#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 13:01:06 2021

@author: andreaapariciomartinez
"""

import numpy as np
import csv
import math as mt
from matplotlib import pyplot as plt
from Network_ID import fn
#%%
chsteps = 70 # parameter goes from 1 to 0 in these many steps

#define simulation data location
folderd = "data/dataSave/"
#nums_tests = np.append(np.array(range(18))+1 ,0)


with open(folderd+fn+'_mean1_o.csv') as csvfile:
                 m1l = list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #mutualistic network
meanA = np.array(m1l) 
with open(folderd+fn+'_samp1_o.csv') as csvfile:
                 s1l = list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #mutualistic network
sampA = np.array(s1l) 

meanAf = np.zeros((len(meanA)+1,meanA.shape[1]))
meanAf[0:len(meanA),:]=meanA

sampAf = np.zeros((len(sampA)+1,sampA.shape[1]))
sampAf[0:len(sampA),:]=sampA

#%%
plt.figure(figsize=(4,8))
plt.subplot(2,1,1)
plt.plot(np.append(np.array(np.arange(1-1/(100*chsteps),-1/(100*chsteps),-1/(100*chsteps))),0),sampAf)
plt.xlim(1,-0.01)
plt.ylabel("abundance")
plt.xlabel("$\mu$")

plt.subplot(2,1,2)
plt.plot(np.append(np.array(np.arange(1-1/chsteps,-1/chsteps,-1/chsteps)),0),meanAf)
plt.xlim(1,-0.01)
plt.ylabel("mean abundance")
plt.xlabel("$\mu$")
#%%
plt.figure(figsize=(6,6))

plt.plot(np.append(np.array(np.arange(1-1/(100*chsteps),-1/(100*chsteps),-1/(100*chsteps))),0),sampAf, alpha=0.1)
plt.ylabel("abundance")
plt.plot(np.append(np.array(np.arange(1-1/chsteps,-1/chsteps,-1/chsteps)),0),meanAf)
plt.xlim(1,-0.01)
plt.xlabel("$\mu$")