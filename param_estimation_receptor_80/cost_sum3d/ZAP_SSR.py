
    
import os, sys, shutil, functools, hashlib, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy.stats
import bionetgen
import threading

from pyswarm import pso
from run_bngl_ZAP import run_bionetgen_ZAP
from interplation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata
from oscar_data_import import expt_import
from ca_bootstrap import bootstrap_ca
from data_slice import data_after_tstart
import sys

def ZAP_only(n,K3,C1,C2,g,time,exp_data,cblb,dir_name):
    
   
    p=2 # which column index which obseravable in bionetgen file #total pZAP file
    
    N=2000 #Ca signal time points and SSR at N timepoints
    
    tstart=1.0 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 25 sec
    
    Vc=25 #pZAP molecules in the simulation box of size 25 um^3
    z=602 #constant factor to convert from molecules/um3 to uM
    
    K2=1
    #print("K1 = "+str(K1)) #SYK concentration
    #print("K2 = "+str(K2)) #Syk proportionality constant

    
    # K1=0.01
    # C1=10000 #Only increases the amplitude 
    # C2=0.1 #As C2 decreases amplitude increases with rising time increase separation of timescale 
    # g=0.05 # the rate of decay : as g increases the amplitude decreases and decay faster
    # m=200
        
    #breakpoint()
#0. Bionetgen parameter set 
    model = bionetgen.bngmodel("BioModel_ZAP70_SYK_v8_test1.bngl")
    model.parameters.Z0 = K3 # setting parameter conc. of ZAP
    model.parameters.SS0 = 0 # setting parameter conc. of SYK
    model.parameters.cblb0 = cblb # setting parameter conc. of SYK
    
    #print(model)

    #print model in directiry name_gamma_HPC.bngl    5 folders
    with open(dir_name+"_Model_ZAP.bngl", "w") as f:
        f.write(str(model)) # writes the changed model to new_model file


   # m.parameters["k"].value = 100 # your new value
   
# #1. Bionetge run : average PZAP   over n run and obseravble is in the (p+1) column
   
    T,avg=run_bionetgen_ZAP(n,p,dir_name) 
    T=np.asarray(T)
    avg=np.asarray(avg)
    
   
    # print(T)
    # print(avg)
    # print(len(T))
    
    
#2.Interpolate PZAP  from 600 points to 2000 points   
    tnew,PZAP=spline(T, avg, N,tstart)
    tnew=np.asarray(tnew)
    PZAP=np.asarray(PZAP) #total number of PZAP molecule in the simulation box of size Vc
    PZAP=(PZAP/(Vc*z)) #pZAP in uM unit
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
    SSR_ZAP=np.sum(np.square(OBJ))
    #print(SSR)
    
    # plot_data=np.column_stack([tnew,ca,idata])
    # file = open('plot_fit.dat', 'w')
    # np.savetxt("plot_fit.dat" , plot_data, fmt=['%f','%f','%f'])
    #plt.plot(TT,MODEL,'g-',tnew,idata,'r-',TT,OBJ,'b--',time,exp_data,'m')
    #plt.show()
    
    
   #8. Delete existing bionetgen output folder and bionetgen model
    for i in range(0,n,1):
        path="/home/gddaslab/ixn004/ZAP_SYK_JAN_2024/FEB_RUNs_2024_kill_CD3z/cost_sum3d/"+dir_name+"_ZAP_HPC_output_"+str(i)
        shutil.rmtree(path)

    os.remove("/home/gddaslab/ixn004/ZAP_SYK_JAN_2024/FEB_RUNs_2024_kill_CD3z/cost_sum3d/"+dir_name+"_Model_ZAP.bngl")
    # for i in range(0,n,1):
    #     path="//Users/ixn004/Dropbox (NCH)/Conference_qbio_Telluride/Cluster/Poster_run/Final_run_cluster_7th_September/Cluster_run/HPC_output_"+str(i)
    #     shutil.rmtree(path)

    #os.remove('/Users/ixn004/Dropbox (NCH)/Conference_qbio_Telluride/Cluster/Poster_run/Final_run_cluster_7th_September/Cluster_run/gamma_HPC.bngl')

   
    #print(time[152])April_parallel_run_23
    #print(exp_data[21])
    print("The (SSR)**0.25 is " + str(SSR_ZAP**(1/4.0)))
    print('remove done')
    
    return SSR_ZAP    