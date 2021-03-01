#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:41:46 2020

@author: andreaapariciomartinez
"""
import csv
import numpy as np
#%%
def buildMat_f(file_n, save=0):

    with open(file_n+'.csv') as csvfile:
             mutnet = list(csv.reader(csvfile,quoting=csv.QUOTE_NONNUMERIC)) #mutualistic network
    MutNet = np.array(mutnet)        
    
    #% build community network from mutualistic matrix

    At = np.zeros((np.sum(MutNet.shape),np.sum(MutNet.shape)))
    
    sp_num = len(At) #number of species
    P_num = MutNet.shape[0] #number of plants (rows)
    A_num = MutNet.shape[1] #number of animals (columns)
    
    # add plant-animal interactions First plants and then animals
    for i in range(P_num):
        for j in range(A_num):
            if MutNet[i,j] != 0:
                At[i,j+P_num] = MutNet[i,j]

    
    # add mutualistic interactions
    
    As= At+At.T #make matrix symmetric
    Abin = np.zeros(As.shape)
    
    #make a binary matrix
    
    for i in range(As.shape[0]):
        for j in range(As.shape[1]):
            if At[i,j]!=0:
                Abin[i,j]=int(1)
            
    
    Abins = Abin+Abin.T
    
    ##%% save files
    if save ==1:
        np.savetxt(file_n+"_Anotsym.csv",Abin)
        np.savetxt(file_n+'_A.csv',As)
        #
        np.savetxt(file_n+'_Ab.csv',Abins)
        
    #% take bidirectionalities away
    #prop = 50
    #bp = mt.floor(len(At)*prop/100)
    #Abp = np.copy(Abin)
    #
    #for i in range(sp_num):
    #    for j in range(sp_num):
    #        if i<j:
    #            if At[i,j] != 0:
    #                if bp>0:
    #                    Abp[j,i]=Abin[i,j]
    #                    bp=bp-1
    #                    print([j,i])
    #        
    #
    #np.savetxt(file_n+'_Ab50.csv',Abp)
    #print(file_n+'_Ab50.csv')
    #

    return [As, Abins]

#%% build competition vector
    
def compMat_f(file_n):
   
    with open(file_n+'.csv') as csvfile:
             mutnet = list(csv.reader(csvfile,quoting=csv.QUOTE_NONNUMERIC)) #mutualistic network
    MutNet = np.array(mutnet)        
    
    P_num = MutNet.shape[0] #number of plants (rows)
    A_num = MutNet.shape[1] #number of animals (columns)
    
    compM = np.zeros((np.sum(MutNet.shape),np.sum(MutNet.shape)))
    
    # add  plant-plant competition
    for i in range(P_num):
        for j in range(A_num):
            if MutNet[i,j] !=0:
                for k in range(i+1,P_num,1):
                    if MutNet[k,j]!=0:
                        compM[i,k]=-1
    
    # add  animal-animal competition
                        
    for j in range(A_num):
        for i in range(P_num):
            if MutNet[i,j]!=0:
                for k in range(j+1,A_num,1):
                    if MutNet[i,k]!=0:
                        compM[j+P_num,k+P_num]=-1    

    return [compM, P_num, A_num]