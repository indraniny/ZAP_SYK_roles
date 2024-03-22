#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 21:33:25 2022

@author: ixn004
"""

#22nd January 2024 ZAP SYK model,  one line corrected at data_slice.py #for SYK only
#8th NOvV -5 parameter estimation for Human
#6 November 2023
#20th May
#7th May, 2023- Mice run remove 3 parameters
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
import threading

from pyswarm import pso
from run_bngl import run_bionetgen
from interplation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata
from oscar_data_import import expt_import
from ca_bootstrap import bootstrap_ca
from data_slice import data_after_tstart
import sys

from run_bngl import run_bionetgen
from ZAP_SSR import ZAP_only
from run_bngl_ZAP import run_bionetgen_ZAP

from SYK_SSR import SYK_only
from run_bngl_SYK import run_bionetgen_SYK




def ca_cost_func(param):
    
    param=10 ** param
    

        #[ZAP0,C1,C2,g,KSA,KST,C1s,C2s,gs]

    K3=param[0] #SYK0 concentration
    C1=param[1] #propensity
    C2=param[2] #uSYK-> pSYK auto
    g=param[3] # pSYK+uSYK -> pSYK+pSYK trans
    KSA = param[4] #estimated using ZAP only
    KST= param[5]
    C1s= param[6]
    C2s= param[7]
    gs= param[8]
    
    n=3

    dir_name=str(outdir)
    #ZAP
    SSR_ZAP=ZAP_only(n,K3,C1,C2,g,time_ZAP,exp_data_ZAP,dir_name)
    
    K1=(K3/1.5)
    #SYK
    SSR_SYK=SYK_only(n,K1,KSA,KST,C1s,C2s,gs,time_SYK,exp_data_SYK,dir_name)



    cost_sum=SSR_ZAP+SSR_SYK
    return cost_sum   
    

def main():
#5. PSO
      
        global outdir 
        outdir=sys.argv[1] #this will contain the string analysis${i} output directories
        #os.mkdir(outdir)

      #2.start by importing mean data from Albert
      
        #Human SYK only   
    
        data_ZAP=pd.read_csv('SYKKO_ca.dat',sep="\t", comment='#', header=None) #ZAP only
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
    
        data_SYK=pd.read_csv('ZAPKO_ca.dat',sep="\t", comment='#', header=None) #SYK only
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
        data_WT=pd.read_csv('WT_ca.dat',sep="\t", comment='#', header=None)
        global time_WT, exp_data_WT
        #if the data has a time column:
        data_WT=np.asarray(data_WT) 
        #print(type(data))
        time_WT=data_WT[:,0]
        #print(len(time))
        #print(time)
        exp_data_WT=data_WT[:,1]
        #plt.plot(time,exp_data,'g--')
        #print(type(time))
        
        

        #[ZAP0,C1,C2,g,KSA,KST,C1s,C2s,gs]
        lb=[2,  -4.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0]
        ub=[3.2, 1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0,  1.0]
        #breakpoint()        
     
        optimized_param, residue=pso(ca_cost_func, lb=lb, ub=ub,maxiter=20, swarmsize=40)
        
        #define a lock aloow one thread at a time, all other thread muct wait until the lock is released


        #optimized_param[2]=10**3.87811026
        #optimized_param[3]=10**-0.27792396 
        #optimized_param[4]=10**-2.03077759

        thlock=threading.Lock()       
        
        print(optimized_param)
        print(residue)
        xx=math.sqrt(residue)
        yy=math.sqrt(xx) 
        
        print(type(optimized_param))
        optimized_param=np.asarray(optimized_param)
      
        f = open(str(outdir)+"_param_residue.dat", 'w')
        f.write("The optimized parameter is "+str(optimized_param)+"\n Residue = "+str(residue)+"\n Squared residue = " +str(xx)+"\n 1/4th residue = "+str(yy))
        f.close()
         
        thlock.acquire() #lock on
        f2 = open("common.dat", 'a') #append in a common file
        f2.write(str(outdir)+"\t"+str(optimized_param)+"\n")
        thlock.release() #lock off
        

                    
        thlock.acquire() #lock on
        f3 = open("common_parameters.dat", 'a')
        f3.write(str(optimized_param[0])+"\t"+str(optimized_param[1])+"\t"+str(optimized_param[2])+"\t"+str(optimized_param[3])+"\t"+str(optimized_param[4])+"\t"+str(optimized_param[5])+"\t"+str(optimized_param[6])+"\t"+str(optimized_param[7])+"\t"+str(optimized_param[8])+"\t"+str(yy)+"\n")
        thlock.release() #lock off
        

if __name__=="__main__" :
    main()    
    
   

    
