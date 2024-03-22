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
      

def ca_cost_func(param):
    
    param=10 ** param
 
    K1=param[0] #SYK0 concentration
    K2=param[1] #propensity
    KSA=param[2] #uSYK-> pSYK auto
    KST=param[3] # pSYK+uSYK -> pSYK+pSYK trans
    C1 = param[4] #estimated using ZAP only
    C2 = param[5]
    g = param[6]

    K1=434

    dir_name=str(outdir)

    # print(dir_name)
    # print(type(dir_name))    
 

    n=1 #average of howmany bionetgen run
    p=1 # which column index which obseravable in bionetgen file
    
    N=2000 #Ca signal time points and SSR at N timepoints
    
    tstart=1.0 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 25 sec
    
    Vc=25 #pZAP molecules in the simulation box of size 25 um^3
    z=602 #constant factor to convert from molecules/um3 to uM
    
    K2=1
    print("K1 = "+str(K1)) #SYK concentration
    print("K2 = "+str(K2)) #Syk proportionality constant

    
    # K1=0.01
    # C1=10000 #Only increases the amplitude 
    # C2=0.1 #As C2 decreases amplitude increases with rising time increase separation of timescale 
    # g=0.05 # the rate of decay : as g increases the amplitude decreases and decay faster
    # m=200
        
    #breakpoint()
#0. Bionetgen parameter set 
    model = bionetgen.bngmodel("BioModel_ZAP70_SYK_v8_test1.bngl")
    model.parameters.SS0 = K1 # setting parameter conc. of SYK
    model.parameters.KSA0 = KSA # setting parameter conc. of SYK
    model.parameters.KST0 = KST # setting parameter conc. of SYK
    
    #print(model)

    #print model in directiry name_gamma_HPC.bngl    5 folders
    with open(dir_name+"_Model_ZAP70_SYK.bngl", "w") as f:
        f.write(str(model)) # writes the changed model to new_model file


   # m.parameters["k"].value = 100 # your new value
   
# #1. Bionetge run : average PZAP   over n run and obseravble is in the (p+1) column
   
    T,avg=run_bionetgen(n,p,dir_name) 
    T=np.asarray(T)
    avg=np.asarray(avg)
    
   
    # print(T)
    # print(avg)
    # print(len(T))
    
    
#2.Interpolate PZAP  from 600 points to 2000 points   
    tnew,PZAP=spline(T, avg, N,tstart)
    tnew=np.asarray(tnew)
    PZAP=np.asarray(PZAP) #total number of PZAP molecule in the simulation box of size Vc
    PZAP=K2*(PZAP/(Vc*z)) #pZAP in uM unit
    # plt.plot(T,avg,'g*',tnew,PZAP,'r-') 
    # plt.show()  
    
  
#3. Take the experimental data and interpolate at the theoretical points
    idata,tnew=exp_idata(time, exp_data, tnew)
    idata=np.asarray(idata)
    tnew=np.asarray(tnew)
    #print(idata)
    
    
#4. Make TT, MODEL, EXPT after tstart      
    TT,EXPT,count, CA0 = data_after_tstart(tnew,idata)
    # print(count)
    # print(CA0)
    # print(len(tnew))
    # print(len(idata))
    
#5. Ca Signal   from the ODE model uisng the PZAP signal   
    y0 = [CA0, 1] #changed
    ca,h=solve(tnew,N,y0,PZAP,C1,C2,g)
    ca=np.asarray(ca)
    h=np.asarray(h)
    # plt.plot(tnew,ca,'m')
    # plt.show()
    
    #4. Print to a file    
    
    #data = np.column_stack([tnew, ca,h,PZAP])
    #datafile_path = "/your/data/output/directory/datafile.txt"
    #np.savetxt("t_ca.dat" , data, fmt=['%f','%f','%f','%f'])
    #file1.close()
   
    
    #4. Plot  
    # plt.plot(tnew,ca, 'b', label='ca(t)')
    # plt.plot(tnew, h, 'g', label='h(t)')
    # plt.plot(tnew, PZAP, 'r:', label='PZAP(t)')
    # plt.legend(loc='best')
    # plt.xlabel('t')
    # plt.grid()
    # plt.show()    
    
    
        
  #6. fit starts from tstart, so, MODEL and EXPT data containts data after tstart    
    MODEL=[]
    for i in range (count,len(tnew),1): 
        MODEL.append(ca[i])
 
    TT=np.asarray(TT)
    EXPT=np.asarray(EXPT)
    MODEL=np.asarray(MODEL)
    #print(TT)
    
  #7. Finding objective function only at t>25 sec
    ca_diff=(idata-ca)
    OBJ=(EXPT-MODEL)
    
    #print(ca_diff)
    SSR=np.sum(np.square(OBJ))
    #print(SSR)
    
    # plot_data=np.column_stack([tnew,ca,idata])
    # file = open('plot_fit.dat', 'w')
    # np.savetxt("plot_fit.dat" , plot_data, fmt=['%f','%f','%f'])
    #plt.plot(TT,MODEL,'g-',tnew,idata,'r-',TT,OBJ,'b--',time,exp_data,'m')
    #plt.show()
    
    
   #8. Delete existing bionetgen output folder and bionetgen model
    for i in range(0,n,1):
        path="/home/gddaslab/ixn004/ZAP_SYK_JAN_2024/FEB_RUNs_2024_kill_CD3z/Cluster_run_SYK_only_konfixed_phosph_C2LCKc/"+dir_name+"_HPC_output_"+str(i)
        shutil.rmtree(path)

    os.remove("/home/gddaslab/ixn004/ZAP_SYK_JAN_2024/FEB_RUNs_2024_kill_CD3z/Cluster_run_SYK_only_konfixed_phosph_C2LCKc/"+dir_name+"_Model_ZAP70_SYK.bngl")
    # for i in range(0,n,1):
    #     path="//Users/ixn004/Dropbox (NCH)/Conference_qbio_Telluride/Cluster/Poster_run/Final_run_cluster_7th_September/Cluster_run/HPC_output_"+str(i)
    #     shutil.rmtree(path)

    #os.remove('/Users/ixn004/Dropbox (NCH)/Conference_qbio_Telluride/Cluster/Poster_run/Final_run_cluster_7th_September/Cluster_run/gamma_HPC.bngl')

   
    #print(time[152])April_parallel_run_23
    #print(exp_data[21])
    print("The (SSR)**0.25 is " + str(SSR**(1/4.0)))
    print('remove done')
    
    return SSR    
    

def main():
#5. PSO
      #K1=0.0005
      #A=ca_cost_func(K1)
      # print(A)
      
      
      
     #1.start by importing your data  and bootstrapping  
        #m=300 # how many bin you want to take from the raw flow file
        #global time, exp_data
        #time,exp_data=bootstrap_ca(m)
      
        global outdir 
        outdir=sys.argv[1] #this will contain the string analysis${i} output directories
        #os.mkdir(outdir)

      #2.start by importing mean data from Albert
      
        #data=pd.read_csv('WT_ca.dat',sep="\t", comment='#', header=None)
        #data=pd.read_csv('ZAPKO_ca.dat',sep="\t", comment='#', header=None)
        data=pd.read_csv('ZAPKO_ca.dat',sep="\t", comment='#', header=None)
        global time, exp_data
        #if the data has a time column:
        data=np.asarray(data) 
        #print(type(data))
        time=data[:,0]
        #print(len(time))
        #print(time)
        exp_data=data[:,1]
        #plt.plot(time,exp_data,'g--')
    
    
        #print(type(time))
        
        
        
        # #parameter bounds [SYK0, C4,Ksa,kst,C1,C2,g]
        # #parameter bounds [SYK0,C4,ksa,kst,C1,C2,g]
    
        #lb=[2.4,-3.0, -3.0,-3.3,-5.0, -5.0, -6.0]
        #ub=[2.9, 0.0, -2.0,-2.3, 1.0, 1.0, 1.0]

        lb=[2,-3.0, -3.0,-3.3,-5.0, -5.0, -5.0]
        ub=[3.0, 0.0, -2.0,-2.3, 1.0, 1.0, 1.0]
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
        f3.write(str(optimized_param[0])+"\t"+str(optimized_param[1])+"\t"+str(optimized_param[2])+"\t"+str(optimized_param[3])+"\t"+str(optimized_param[4])+"\t"+str(optimized_param[5])+"\t"+str(optimized_param[6])+"\t"+str(yy)+"\n")
        thlock.release() #lock off
        

if __name__=="__main__" :
    main()    
    
   

    
