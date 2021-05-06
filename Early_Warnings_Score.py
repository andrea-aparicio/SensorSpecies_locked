#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 11:19:28 2020

@author: andreaapariciomartinez
"""
import scipy.stats as sts
import numpy as np
import csv
import math as mt
from matplotlib import pyplot as plt
from Network_ID import fn
#%%
chsteps = 70 # parameter goes from 1 to 0 in these many steps

#define simulation data location
folderd = "data/dataSave/"
nums_tests = [13,14,15,19]
for k in nums_tests:

#get simulation data
    with open(folderd+fn+'_var1_'+str(k)+'.csv') as csvfile:
                 v1l = list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #mutualistic network
    varM100 = np.array(v1l) 
    sp = varM100.shape[1]
#    with open(folderd+fn+'_var0.csv') as csvfile:
#                 v0l = list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #mutualistic network
#    varM0 = np.array(v0l)   
#    sp = len(varM0.T)    
#    chsteps = len(varM0)
    
#% calculate slope for a rolling window
    win = mt.ceil(chsteps/3)
#    polyf0 = np.zeros((chsteps-win,sp))
    polyf100 = np.zeros((chsteps-win,sp))
    a=range(win)
    for i in range(chsteps-win):
        for j in range(sp):
#            polyf0[i,j] = np.polyfit(a,varM0[i:i+win,j],1)[0]
            polyf100[i,j] = np.polyfit(a,varM100[i:i+win,j],1)[0]
  
#% calculate Early Warning Score (when slope is never negative again)
#    detection0t = win*np.ones(sp)
    detection100t = win*np.ones(sp)


#    for j in range(sp):
#        for i in range(chsteps-win):
#            if polyf0[chsteps-win-1-i,j]<=0:
#                detection0t[j] = chsteps-i
#                break
    for j in range(sp):
        for i in range(chsteps-win):
            if polyf100[chsteps-win-1-i,j]<=0:
                detection100t[j] = chsteps-i
                break
#    detection0 = 1-detection0t/chsteps
    detection100 = 1-detection100t/chsteps
    plt.scatter(range (sp),detection100)   
#    np.savetxt(folderd+fn+"_detection0.csv", detection0)
    np.savetxt(folderd+fn+'_detection100_'+str(k)+'.csv', detection100)
    
    
