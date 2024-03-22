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
    

    
# =============================================================================
    
    K0=0.0005
    
  
    K1=0 #SYK concentration
    K2=0
    K3=10**(2.911034846) #ZAP concentration
    
    KSA=10**(-3.466526711)
    KST=10**(-3.955852847)
    
    C1 = 10**(-0.64877186)
    C2 = 10**(-2.194741228)
    g =10**(-1.8403492)



    n=1#average of howmany bionetgen run
    p=1 # which column index which obseravable in bionetgen file SYK
    p1=4 #ZAP
    
    N=2500 #Ca signal time points and SSR at N timepoints
    
    tstart=1.0 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 60 sec
    
    Vc=25 #pZAP molecules in the simulation box of size 25 um^3
    z=602 #constant factor to convert from molecules/um3 to uM
    
    print("K1 = "+str(K1))
    print("K2 = "+str(K2))
    print("K3 = "+str(K3))
    print("C1 = "+str(C1))
    print("C2 = "+str(C2))
    print("g = "+str(g))


        
    #breakpoint()
#0. Bionetgen parameter set 

    #model = bionetgen.bngmodel("BioModel_ZAP70_SYK_v8_test1d.bngl") #homodimer-gamma code
    model = bionetgen.bngmodel("BioModel_ZAP70_SYK_v8_test1e.bngl") #homodimer-gamma code
    
    
    model.parameters.SS0 = K1 # setting parameter conc. of SYK
    model.parameters.Z0 = K3 # setting parameter conc. of ZAP
    model.parameters.KSA0 = KSA # setting parameter conc. of SYK
    model.parameters.KST0 = KST # setting parameter conc. of SYK
  
	    
    
    #print(model)

    #print model in directiry name_gamma_HPC.bngl    5 folders
    with open("Model_ZAP70_SYK.bngl", "w") as f:
        f.write(str(model)) # writes the changed model to new_model file


   # m.parameters["k"].value = 100 # your new value
   
# #1. Bionetge run : average PZAP   over n run and obseravble is in the (p+1) column
   
    T,avg_SYK=run_bionetgen(n,p) #SYK
    T=np.asarray(T)
    avg_SYK=np.asarray(avg_SYK)#SYK signal
    
 
    # print(T)
    # print(avg)
    # print(len(T))

    
#2.Interpolate pSYKP  from 600 points to 2000 points   
    tnew,PSYK=spline(T, avg_SYK, N,tstart)
    tnew=np.asarray(tnew)
    PSYK=np.asarray(PSYK) #total number of PZAP molecule in the simulation box of size Vc
    plt.plot(T,avg_SYK,'g*',tnew,PSYK,'r-') 
    
    
        
    T,avg_ZAP=run_bionetgen(n,p1) #SYK
    T=np.asarray(T)
    avg_ZAP=np.asarray(avg_ZAP)#ZAP signal
    
        
#2.Interpolate PZAP  from 600 points to 2000 points   
    tnew,PZAP=spline(T, avg_ZAP, N,tstart)
    tnew=np.asarray(tnew)
    PZAP=np.asarray(PZAP) #total number of PZAP molecule in the simulation box of size Vc
    plt.plot(T,avg_ZAP,'m*',tnew,PZAP,'g-') 
    
    PZAP=PZAP/(Vc*z) #pZAP in uM unit
    PSYK=PSYK/(Vc*z)

    plt.show()  
    
    #PZAP interpolation on avg is very good 22nd Jan 2024    
    
  
#3. Take the experimental data and interpolate at the theoretical points
    idata,tnew=exp_idata(time, exp_data, tnew)
    idata=np.asarray(idata)
    tnew=np.asarray(tnew)
    #print(idata)
    
    
#4. Make TT, MODEL, EXPT after tstart      
    TT,EXPT,count, CA0 = data_after_tstart(tnew,idata)
    print(f'count=',count)
    #print(CA0)
    #print(EXPT)
    # print(len(tnew))
    # print(len(idata))
    
# # #5. Ca Signal from the ODE model uisng the PZAP signal   ##########PZAP
   
    
    y0z = [CA0, 1]
    ca_PZAP,h_PZAP=solve(tnew,N,y0z,PZAP,C1,C2,g)
    ca_PZAP=np.asarray(ca_PZAP)
    h_PZAP=np.asarray(h_PZAP)
    #plt.plot(tnew,ca,'m') #model generated
    print(f'ca0 in main',ca_PZAP[0])
    # plt.show()
    print(ca_PZAP)
    

    

    
  
    ca=ca_PZAP

    
  #6. fit starts from tstart, so, MODEL and EXPT data containts data after tstart    
    MODEL=[]
    for i in range (count,len(tnew),1): 
        MODEL.append(ca[i])
 
    TT=np.asarray(TT)
    EXPT=np.asarray(EXPT)
    MODEL=np.asarray(MODEL)
    #print(TT)
    #plt.plot(TT,EXPT,'b*',time,exp_data,'r-')
    #plt.plot(TT,MODEL, 'g*',TT,EXPT,'b-',time,exp_data,'r-')
    
    #TT, EXPT are the interpolated data for time > 60 sec
    
  #7. Finding objective function only at t>60 sec
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
    
  
    
    plot_data = np.column_stack([TT, MODEL, EXPT])
    file = open('plot_fit.dat', 'w')
    
    #np.savetxt("plot_fit.dat", plot_data, fmt=['%f', '%f', '%f'])
    
    #np.savetxt("mgamma_model_ca_prediction.dat", plot_data, fmt=['%f', '%f', '%f'])
    #np.savetxt("hCD3z_model_ca_prediction.dat", plot_data, fmt=['%f', '%f', '%f'])
    #np.savetxt("WTz_model_ca_prediction.dat", plot_data, fmt=['%f', '%f', '%f'])
    
    plt.rc('xtick', labelsize=30)
    plt.rc('ytick', labelsize=30)
    # plt.ylim(20000,75000)
    # plt.plot(TT, MODEL, 'c-', tnew, idata, 'r-',
       #     TT, OBJ, 'b--', time, exp_data, 'm')
    #plt.plot(tnew,idata,'c',TT,MODEL,'k',linewidth=3.0)
    

    plt.plot(tnew,ca,'c-',TT,MODEL, 'g-',TT,EXPT,'b*',time,exp_data,'r-',linewidth=2.0)
    #plt.plot(TT,MODEL,'k',linewidth=3.0)
    #plt.show()

    #cmap = plt.cm.Reds
   
    
    # Set the y-axis ticks between 0.2 and 1 with a gap of 0.1
    y_ticks = np.arange(0.2, 1.0, 0.1)
    plt.yticks(y_ticks)
   # plt.plot(TT,EXPT,color=cmap(0.8),linewidth=3.0)
    plt.plot(TT,EXPT,'skyblue',linewidth=3.0)
  
    #plt.plot(TT,EXPT,'blue',linewidth=3.0)
    
    
    plt.plot(TT,MODEL,'k',linewidth=3.0)
    #plt.show()
    # plt.plot(TT,EXPT,"k-",linewidth=3.0)
    #plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    # Create a boxplot and set the width of the boxes
    # Customize the axes line width
    ax = plt.gca()  # Get the current axes
    ax.spines['left'].set_linewidth(2.0)  # Adjust the left spine's line width
    ax.spines['bottom'].set_linewidth(2.0)  # Adjust the bottom spine's line width
    ax.spines['right'].set_linewidth(2.0)  # Adjust the right spine's line width
    ax.spines['top'].set_linewidth(2.0)  # Adjust the top spine's line width
    plt.show()
    
    
    print("The (SSR)**0.25 is SSR=" + str(SSR**(1/4.0)))
    
    print("K1 = "+str(K1))
    print("K2 = "+str(K2))
    print("K3 = "+str(K3))
    print("KSA = "+str(KSA))
    print("KST = "+str(KST))
    print("C1 = "+str(C1))
    print("C2 = "+str(C2))
    print("g = "+str(g))
    print("SSR=" + str(SSR**(1/4.0)))
    
    
    print('remove done')

    
    return SSR    
    

def main():
#5. PSO
     
        
        #Human List
        #data=pd.read_csv('humanCD16_hCD3z.dat',sep="\t", comment='#', header=None)
        #data=pd.read_csv('ZAPKO_ca.dat',sep="\t", comment='#', header=None) #SYK only
        data=pd.read_csv('SYKKO_ca.dat',sep="\t", comment='#', header=None) #ZAP only
        #data=pd.read_csv('WT_ca.dat',sep="\t", comment='#', header=None)
        
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
        
        K0=0.0005
        A=ca_cost_func(K0)
        

if __name__=="__main__" :
    main()    
    
   

    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 12:32:52 2024

@author: ixn004
"""

