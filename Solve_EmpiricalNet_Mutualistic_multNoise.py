#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 17:09:28 2020

@author: andreaapariciomartinez
"""
import csv
import numpy as np
import sdeint
from matplotlib import pyplot as plt
import pandas as pd
from Network_ID import fn
#%% time and simulation parameters
chsteps = 10 # parameter goes from 1 to 0 in these many steps
#-P..chsteps = 70 # parameter goes from 1 to 0 in these many steps
#-P..tch = 50 # sims run for this amount of time before changing the parameter
tch = 20 # sims run for this amount of time before changing the parameter
tini = 0.0 #initial time
st = .01 #step size
tmax= tch*chsteps+(2*tch) #final time 
tarray=np.arange(tini,tmax,st)
dch = np.zeros(len(tarray)) #stores the values of the parameter for plots
for i in range(len(dch)):
    dch[i] = 1-((1/chsteps)*np.floor((tarray[i])/tch))
dchc=np.arange(1,0,-1/len(tarray))
#-P..p_start=1 #value to start the decreasing parameter
p_start=.5 #value to start the decreasing parameter
selfW = -1
h = .01 #handling time
#P-..sigmax = 10 #additive noise gain 
sigmax = .5 # noise gain 
#%% Function to change a parameter from p_start to 0, in chsteps of tch time length
def calc_d(t):
    d = p_start-((p_start/chsteps)*np.floor((t)/tch)) 
    if d<0:
        d=0
    return d

#%% Mutualistic dynamics
def Mut_eq(x,t):    
    dx = np.zeros(len(A))
   
    d = calc_d(t) #parameter that changes

    x[x<0]=0
    
    for i in range(len(A)): #system dynamics
        dx[i] = x[i]*(r[i]+selfW*x[i]-Ac[i,:]@x+((d*A[i,:]@x)/(1+h*d*A[i,:]@x)))#- 10e-4
    
    return dx
    
#%% multiplicative noise 
def mutN_eq(x,t):
    gv=sigmax*np.ones(len(x))*x
    return np.diag(gv)

#%% define data destination folder    
folder = "Community_parameters/"  #parameter folder
folderd = "data/" #data destination folder

#%%  build system matrix from parameters files 
nums_tests = [0]
for i in nums_tests:
    print(["test",i])
    num_test = str(i)
    with open(folder+fn+'_A_'+num_test+'.csv') as csvfile:
                 Al = list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #mutualistic network
    A = np.array(Al)      

    with open(folder+fn+'_Ac_'+num_test+'.csv',newline='') as csvfile:
             Acl =  list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #competition network
    Ac = np.array(Acl)

    with open(folder+fn+'_r_1.csv',newline='') as csvfile:
                 rl =  list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #growth rate
    r = np.array(rl).flatten()

    with open(folder+fn+'_x0_1.csv',newline='') as csvfile:
                 x0l =  list(csv.reader(csvfile,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)) #initial conditions
    x0 = np.array(x0l).flatten()
    
    print(fn)


##%% solve stochastic system
    Muts1 = sdeint.itoSRI2(Mut_eq, mutN_eq, x0, tarray)
    print("100done")
#    Ac=np.zeros(A.shape) #remove competition
#    Muts0 = sdeint.itoSRI2(Mut_eq, addN_eq, x0, tarray)
#    print("0done")
    
##%% plot to check simulation results
    plt.figure()
    plt.plot(Muts1)
    plt.ylim(0,1)
    plt.savefig(fn+"_ts1_"+num_test+".png", dpi=300)
#    plt.close()
    
#    plt.figure()
#    plt.plot(Muts0)
#    plt.ylim(0,100)
#    plt.savefig(fn+"_ts0.png", dpi=300)
#    plt.close()
    
##%% calculate variance, means and take a sample of timeseries
#    posr=100
#    pextra = 200 #skips these steps after every change to avoid transitory jumps
#    int(tch/st)-1
#    ind = np.arange(0,int(tch/st)-1-pextra,int((int(tch/st)-1-pextra)/posr+1))
#    
#    varM = np.zeros((int(chsteps),len(Muts1.T)))
#    meanM = np.zeros((int(chsteps),len(Muts1.T)))
#    sampM = np.zeros((int(chsteps*posr),len(Muts1.T)))
#    
##    varM0 = np.zeros((int(chsteps),len(Muts0.T)))
##    meanM0 = np.zeros((int(chsteps),len(Muts0.T)))
##    sampM0 = np.zeros((int(chsteps*posr),len(Muts0.T)))
#    
#    for i in range(int(chsteps)):
#        for j in range(len(Muts1.T)):
#            tempM = pd.Series(Muts1[i*int(tch/st)+pextra:(i*int(tch/st))+int(tch/st)-1,j])
#            varM[i,j] = pd.DataFrame.var(tempM)
#            meanM[i,j] = pd.DataFrame.mean(tempM)
#            sampM[(i*posr):(i*posr)+posr,j] =  np.array(tempM[ind])
#            
#            tempM0 = pd.Series(Muts0[i*int(tch/st)+pextra:(i*int(tch/st))+int(tch/st)-1,j])
#            varM0[i,j] = pd.DataFrame.var(tempM0)
#            meanM0[i,j] = pd.DataFrame.mean(tempM0)
#            sampM0[(i*posr):(i*posr)+posr,j] =  np.array(tempM0[ind])
#            

#% save data
    
#    np.savetxt(folderd+fn+"_var1_"+num_test+".csv", varM)
#    np.savetxt(folderd+fn+"_mean1_"+num_test+".csv", meanM)
#    np.savetxt(folderd+fn+"_samp1_"+num_test+".csv", sampM)
#    print(["saved",i])
##    
#    np.savetxt(folderd+fn+"_var0.csv", varM0)
#    np.savetxt(folderd+fn+"_mean0.csv", meanM0)
#    np.savetxt(folderd+fn+"_samp0.csv", sampM0)
#    
