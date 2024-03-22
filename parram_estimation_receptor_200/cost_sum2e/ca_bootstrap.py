#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 12:08:40 2022

@author: ixn004
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import FlowCal #import module

import random

def bootstrap_ca(n):

    A=pd.read_csv('raw_ca.csv')
    # print(A.shape)
    
    A=np.asarray(A)
    # print(A.shape)
    
    ca=A[:,13]
    t=A[:,14]
    # plt.plot(t,ca)
    # plt.show()
    
    #plt.hist2d(t,ca,100,cmap='jet')
    
    # t=np.asarray(t)
    # ca=np.asarray(ca)
    
    t=list(t)
    ca=list(ca)
  
    
    #n=100.0
    tmax=max(t)
    dt=tmax/n #data is divided into the number of bins
    
    
    length=len(t)
    #print(length)
    
    t1=0.0
    t2=t1+dt
    k=0
    j=0
    mean_t=[] #define a list
    mean_ca=[] #define a list
    T=[]
    j=0
    X=[] #define a list
    Y=[] #define a list
    ca_boot=[]
    
    Y1=[]
    Y2=[]
    
    
    my_samples = []
    
    for i in range(length):
       
        
        if t[i] > t1 and t[i] < t2:
         
            X.append(t[i])
            Y.append(ca[i])
            #print(Y,len(Y))
            
    
   
    #print(Y)
        elif t[i] > t2: 
            #print(str(X)+'done')
            
            #YY=random.sample(Y, 200) ##create a sample of size 200
            
            YY = np.random.choice(Y, size=len(Y), replace=True)
            my_samples.append(YY.mean())
            
            #print(len((X)))
            X=np.asarray(X)
            mean_t.append(X.mean())
            
            
            #print(mean_t)
            # mean_ca.append(sum(Y)/len(Y))
            # ca_boot.append(sum(YY)/len(YY))
            #print(mean_t,mean_ca)
            
            X=[] #define a list
            Y=[] #define a list
            YY=[]
            t1=t2
            t2=t1+dt
        
     
    
        
    # plt.plot(mean_t,ca_boot,'r')
    # plt.show()
    plt.plot(mean_t,my_samples,'g')
    plt.show()
 
    
    return mean_t,my_samples
 # data=pd.read_csv('oscar_ca.dat',sep="\t", comment='#', header=None)
      # #if the data has a time column:
      # data=np.asarray(data)
      # #print(type(data))
      # global time, exp_data
      
      # time=data[:,0]
      # #print(len(time))
      # #print(time)
      # exp_data=data[:,1]  
