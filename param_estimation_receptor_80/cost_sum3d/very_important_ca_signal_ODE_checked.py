#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 21:27:35 2022

@author: ixn004
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 21:33:25 2022

@author: ixn004
"""
#Very important matched with ode23 bay matlab considering PZAP as a unit vector
#All the functions checked on 5th September
# Ca ODE is checked on 6th September

import os, sys, shutil, functools, hashlib, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy.stats
import bionetgen

from pyswarm import pso
from run_bngl import run_bionetgen
from interplation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata
from oscar_data_import import expt_import
from ca_bootstrap import bootstrap_ca
      
def main():
    
   
    n=1 #average of howmany bionetgen run
    p=1 # which column index which obseravable in bionetgen file
    
    N=2000 #Ca signal time points and SSR at N timepoints
    
    tstart=1 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 25 sec
    
    Vc=25 #pZAP molecules in the simulation box
    z=602 #constant factor to convert from molecules/um3 to uM
    
    K1=0.01
    C1=10000 #Only increases the amplitude 
    C2=0.1 #As C2 decreases amplitude increases with rising time increase separation of timescale 
    g=0.05
    m=200
    
    
    global time, exp_data
    time,exp_data=bootstrap_ca(m)
    

    
#0. Bionetgen parameter set 
    model = bionetgen.bngmodel("gamma_HPC_0.bngl")
    model.parameters.Kab = K1 # setting parameter k to 1
    #print(model)
    with open("gamma_HPC.bngl", "w") as f:
        f.write(str(model)) # writes the changed model to new_model file


   # m.parameters["k"].value = 100 # your new value
   
# #1. Bionetge run : average PZAP   
   
    T,avg=run_bionetgen(n,p) 
    T=np.asarray(T)
    avg=np.asarray(avg)
    
   
    # print(T)
    # print(avg)
    print(len(T))
    
    
#2.Interpolate PZAP  from 600 points to 2000 points   
    tnew,PZAP=spline(T, avg, N,tstart)
    tnew=np.asarray(tnew)
    PZAP=np.asarray(PZAP) #total number of PZAP molecule in the simulation box of size Vc
    PZAP=PZAP/(Vc*z) #pZAP in uM unit
    # plt.plot(T,avg,'g*',tnew,PZAP,'r-') 
    # plt.show()  
  
    
#3. Ca Signal   from the ODE model    
    y0 = [40000, 10]
    ca,h=solve(tnew,N,y0,PZAP,C1,C2,g)
    ca=np.asarray(ca)
    h=np.asarray(h)
    # plt.plot(tnew,ca,'m')
    # plt.show()
    
  
    
#4. Print to a file    
    
    data = np.column_stack([tnew, ca,h,PZAP])
    #datafile_path = "/your/data/output/directory/datafile.txt"
    #np.savetxt("t_ca.dat" , data, fmt=['%f','%f','%f','%f'])
    #file1.close()
    #4. Plot  
    plt.plot(tnew,ca, 'b', label='ca(t)')
    #plt.plot(tnew, h, 'g', label='h(t)')
    plt.plot(tnew, PZAP, 'r:', label='PZAP(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()
    
    #5. Take the experimental data and interpolate at the theoretical points
    idata,tnew=exp_idata(time, exp_data, tnew)
    idata=np.asarray(idata)
    tnew=np.asarray(tnew)
    #print(idata)
    
          
       
   

    ca_diff=(idata-ca)
    #print(ca_diff)
    SSR=np.sum(np.square(ca_diff))
    #print(SSR)
    
    plot_data=np.column_stack([tnew,ca,idata])
    file = open('plot_fit.dat', 'w')
    np.savetxt("plot_fit.dat" , plot_data, fmt=['%f','%f','%f'])
    plt.plot(tnew,ca,'g-',tnew,idata,'r-',tnew,ca_diff,'b--')
    plt.show()
    
    

    

if __name__=="__main__" :
    main()    
    
   