#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 20:27:48 2020

@author: andreaapariciomartinez
"""

import numpy as np
import scipy.integrate as scint
from matplotlib import pyplot as plt
import math as mt
from build_community_matrix import buildMat_f
from build_community_matrix import compMat_f
from Network_ID import fn
#%% time and simulation parameters  and equations to check realization
tini = 0.0 #initial time
st = .1#.01 #step size
tmax=50
tarray=np.arange(tini,tmax,st)

#% Mutualistic dynamics
    
def Mut_eq(x,t):    
    dx = np.zeros(len(A))
   
    d=1
    x[x<0]=0
    
    for i in range(len(A)): #system dynamics
        dx[i] = x[i]*(r[i]+selfW*x[i]-Ac[i,:]@x+((d*A[i,:]@x)/(1+h*d*A[i,:]@x)))#- 10e-4
    
    return dx

#%% system parameters
h = .01 #handling time
tau = 15
g0 = 8*tau #average level of mutualistic strenght
delta = .5 #trade-off modulation


    folder="WoL_matrices/" #path to folder that stores WoL csv network matrix

    check = 1
    while check!=0:
        print(filename)
        name = folder+filename
        [Aval, Abin] = buildMat_f(name,0) # get interaction matrix
        sp_n = len(Abin)
        r_sp = range(sp_n)  
        [compM, pn, an] = compMat_f(name) #get competition matrix and number of plants and animals
    ##%% Add weights to interactions matrix
        At = np.zeros(Abin.shape)
        for i in r_sp:
            for j in r_sp:
                if Abin[i,j]!=0:
                    At[i,j] = g0*Abin[i,j]/((sum(Abin[i,:])/2)**delta)
    ##%% Add weights to competition matrix
        Ac = np.zeros(compM.shape)
        Bp = 1/pn
        Ba = 1/an
        for i in range(pn):
            for j in range(pn):
                if compM[i,j]!=0:
                    Ac[i,j] =  np.random.uniform(.005, 2*Bp)
        for i in range(pn,sp_n,1):
            for j in range(pn,sp_n,1):
                if compM[i,j]!=0:
                    Ac[i,j] =  np.random.uniform(.005, 2*Ba) 

    ##%% get random initial conditions 
        x0m = 0 #minimum
        x0M = 10 #maximum
        x0 = (x0M-x0m)*np.random.random_sample(len(At))+x0m
        x0 = x0.T #transposing to make x0 is a column array
    ##%% get random growth rate
        rmin = -.5
        rmax = -.1
        r = (rmin-rmax)*np.random.random_sample(len(At))+rmin 
    
    ##%% Check realization
        selfW=-1
        A=At #determined weights (as in Dakos PNAS14)
    ##%% solve deterministic system 
        Mut = scint.odeint(Mut_eq, x0, tarray,  hmax=st) 
    ##%% plot to check simulation results
    
        plt.plot(Mut)
        plt.ylim(-1,100)

        check = sum(Mut[-1,:]==0)
        print(check)
##%% save parameters
        
        folderd = 'Community_parameters/'
        np.savetxt(folderd+filename+"_x0.csv",x0)
        np.savetxt(folderd+filename+"_r.csv",r)
        np.savetxt(folderd+filename+"_Ac.csv",Ac)
        np.savetxt(folderd+filename+"_A.csv",At)
        np.savetxt(folderd+filename+"_Abin.csv",Abin)


