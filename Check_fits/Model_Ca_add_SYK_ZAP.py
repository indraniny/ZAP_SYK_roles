#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 21:33:25 2022

@author: ixn004
"""

#8th November Plot for paper, one line corrected at data_slice.py #estimating only SYK concentration concentration
#20th Sept 2023 NK meeting plot
#11th May 2023- run mean
#8th May, 2023- Mice run remove 3 parameters
#20th April, 2023 - Koff [0,1], bootstraap
#19th April. 2023 giving run for constant c1, c2,g
#14th September - thread lock is used 
#14th september - code is ready for parallel run
#All the functions checked on 5th September
# Ca ODE is checked on 6th September, where if I make PZAP signal is a unit vector, it matches with the the PZAP signal using MATLAB ode23
# here is the link to include link for CaODE in python https://apmonitor.com/pdc/index.php/Main/SolveDifferentialEquations
#7th September the objective function is considered after 25 sec



import os, sys, shutil, functools, hashlib, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy.stats
import bionetgen

from pyswarm import pso
from interplation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata

from data_slice import data_after_tstart


from ZAP_SSR import ZAP_only
from run_bngl_ZAP import run_bionetgen_ZAP

from SYK_SSR import SYK_only
from run_bngl_SYK import run_bionetgen_SYK

from WT_SSR import WT
from run_bngl_WT import run_bionetgen_WT

def ca_cost_func():
    
    #param=10 ** param
    
# =============================================================================
  
    
    n=5#average of howmany bionetgen run




######  ZAP only parameters
    
    # K3=10**(2.95264) #ZAP concentration
    # C1 = 10**(-0.649104)
    # C2 = 10**(-0.268880)
    # g =10**(-2.251950)
   
   
   
    # K1=(K3/1.5)
    # KSA=10**(-1.8250723)
    # KST=10**(-4.5491501)
   

   
    # C1s=10**(-1.262904)
    # C2s=10**(-2.97084)
    # gs=10**(-3.5771803)

   
    K3=10**(3.053297371) #ZAP concentration
    C1 = 10**(-0.104816081)
    C2 = 10**(-0.791708367)
    g =10**(-1.924856019)
    
    
    
    K1=(K3/1.5)
    KSA=10**(-2.578147816)
    KST=10**(-2.662747556)
    

    
    C1s=10**(-1.235574311)
    C2s=10**(-3.965926482)
    gs=10**(-2.987637772)


    
    #ZAP
    SSR_ZAP,tnew_ZAP,ca_ZAP,TT_ZAP,MODEL_ZAP,TT_ZAP,EXPT_ZAP=ZAP_only(n,K3,C1,C2,g,time_ZAP,exp_data_ZAP)
    
    
    
    # #SYK
    SSR_SYK,tnew_SYK,ca_SYK,TT_SYK,MODEL_SYK,TT_SYK,EXPT_SYK=SYK_only(n,K1,KSA,KST,C1s,C2s,gs,time_SYK,exp_data_SYK)

    # #cost_sum=(SSR_SYK)+(SSR_ZAP)+SSR_WT
      

    
    # #WT
    SSR_WT,tnew_WT,ca_WT,TT_WT,MODEL_WT,TT_WT,EXPT_WT=WT(n,K3,C1,C2,g,K1,KSA,KST,C1s,C2s,gs,time_WT,exp_data_WT)
    

    #cost_sum=SSR_ZAP
    cost_sum=(SSR_SYK)+(SSR_ZAP)+SSR_WT
    
    #cost_sum=SSR_ZAP
    #print(f'SSR_ZAP and SSR_SYK',SSR_ZAP**0.25,SSR_SYK**0.25)

    return cost_sum
    
        
    

def main():
#5. PSO
      #K1=0.0005
      #A=ca_cost_func(K1)
      # print(A)
    
        #Human ZAP only
        
        data_ZAP=pd.read_csv('SYKKO_ca1.dat',sep="\t", comment='#', header=None) #ZAP only
        
        global time_ZAP, exp_data_ZAP
        #if the data has a time column:
        data_ZAP=np.asarray(data_ZAP) 
        #print(type(data))
        time_ZAP=data_ZAP[:,0]
        #print(len(time))
        #print(time)
        exp_data_ZAP=data_ZAP[:,1]
        #plt.plot(time,exp_data,'g--')
        
        #Human SYK only   
    
        data_SYK=pd.read_csv('ZAPKO_ca1.dat',sep="\t", comment='#', header=None) #SYK only
        global time_SYK, exp_data_SYK
        #if the data has a time column:
        data_SYK=np.asarray(data_SYK) 
        #print(type(data))
        time_SYK=data_SYK[:,0]
        #print(len(time))
        #print(time)
        exp_data_SYK=data_SYK[:,1]
        #plt.plot(time,exp_data,'g--')
        
        
        #Human WT only   
        data_WT=pd.read_csv('WT_ca1.dat',sep="\t", comment='#', header=None)
        global time_WT, exp_data_WT
        #if the data has a time column:
        data_WT=np.asarray(data_WT) 
        #print(type(data))
        time_WT=data_WT[:,0]
        #print(len(time))
        #print(time)
        exp_data_WT=data_WT[:,1]
        #plt.plot(time,exp_data,'g--')
        
   
        ca_cost_func()
        

if __name__=="__main__" :
    main()    
    
   

    
