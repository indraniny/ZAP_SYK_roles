#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 21:32:03 2024

@author: ixn004
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 20:48:14 2024

@author: ixn004
"""


import os, sys, shutil, functools, hashlib, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy.stats
import bionetgen
import threading

from pyswarm import pso
from run_bngl_SYK import run_bionetgen_SYK
from interplation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata

from data_slice import data_after_tstart
import sys


def SYK_only(n,K1,KSA,KST,C1s,C2s,gs,time,exp_data):
    
    
    
    #ZAP param
    K2=0
    K3=0
    C1 = 0
    C2 = 0
    g =0
    
    
    oo=16

    p=1 # which column index which obseravable in bionetgen file SYK
    p1=4 #ZAP
    
    N=2500 #Ca signal time points and SSR at N timepoints
    
    tstart=1.0 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 60 sec
    
    Vc=25 #pZAP molecules in the simulation box of size 25 um^3
    z=602 #constant factor to convert from molecules/um3 to uM
    
    print("K1 = "+str(K1))
    print("K2 = "+str(K2))
    print("K3 = "+str(K3))
  
        
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
    with open("Model_SYK.bngl", "w") as f:
        f.write(str(model)) # writes the changed model to new_model file


   # m.parameters["k"].value = 100 # your new value
   
# #1. Bionetge run : average PZAP   over n run and obseravble is in the (p+1) column
   
    T,avg_SYK=run_bionetgen_SYK(n,p) #SYK
    T=np.asarray(T)
    avg_SYK=np.asarray(avg_SYK)#SYK signal
    
 
    # print(T)
    # print(avg)
    # print(len(T))

    
#2.Interpolate pSYKP  from 600 points to 2000 points   
    tnew,PSYK=spline(T, avg_SYK, N,tstart)
    tnew=np.asarray(tnew)
    PSYK=np.asarray(PSYK) #total number of PZAP molecule in the simulation box of size Vc
    #plt.plot(T,avg_SYK,'g*',tnew,PSYK,'r-') 
    
    
        
    T,avg_ZAP=run_bionetgen_SYK(n,p1) #SYK
    T=np.asarray(T)
    avg_ZAP=np.asarray(avg_ZAP)#ZAP signal
    
        
#2.Interpolate PZAP  from 600 points to 2000 points   
    tnew,PZAP=spline(T, avg_ZAP, N,tstart)
    tnew=np.asarray(tnew)
    PZAP=np.asarray(PZAP) #total number of PZAP molecule in the simulation box of size Vc
    #plt.plot(T,avg_ZAP,'m*',tnew,PZAP,'g-') 
    
    
    
    
    # Create the folder if it doesn't exist
    folder_name = "paper_plot"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # Write the data to a file in the folder
    file_path = os.path.join(folder_name, "SYK_only.txt")
    with open(file_path, 'w') as file:
        for i in range(len(tnew)):
            file.write(f"{tnew[i]} {PZAP[i]} {PSYK[i]}\n")

    print("Data has been written to:", file_path)
    
    
    
    
        
    #Plot ZAP SYK number
    y_ticks = np.arange(0, 3000, 1000)
    plt.yticks(y_ticks)
    plt.plot(tnew, PZAP, color='red', linewidth=6, label='pZAP')   
    plt.plot(tnew, PSYK, color='blue', linewidth=6, label='pSYK')
    #plt.plot(tnew,PSYK,'b',linestyle='--', linewidth=4,label='pSYK') 
    plt.tick_params(axis='both', which='major', labelsize=oo)  # Adjust label size
    plt.tick_params(axis='both', which='major', width=2, length=8)
    ax = plt.gca()  # Get the current axes
    ax.spines['left'].set_linewidth(2.0)  # Adjust the left spine's line width
    ax.spines['bottom'].set_linewidth(2.0)  # Adjust the bottom spine's line width
    ax.spines['right'].set_linewidth(2.0)  # Adjust the right spine's line width
    ax.spines['top'].set_linewidth(2.0)  # Adjust the top spine's line width
    #plt.legend()
    plt.show()

    
    
    
    PZAP=PZAP/(Vc*z) #pZAP in uM unit
    PSYK=PSYK/(Vc*z)
    # #plt.plot(T,avg/(Vc*z),'g*',tnew,PZAP,'r-') 
    # plt.show()  
    
    #PZAP interpolation on avg is very good 22nd Jan 2024

    
    
    #Defining PSYK as PZAP or total signal
    #TOTAL=((PZAP/(Vc*z))+((K2*(PSYK/(Vc*z)) )))/(1+K2)#total signal in uM unit
    #TOTAL=(PZAP+K2*PSYK)/(Vc*z)#total signal in uM unit
    #plt.plot(T,avg/(Vc*z),'g*',tnew,PZAP,'r-') 

    
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

    
    CA0z=0.24724    
    y0z = [CA0, 1]
    ca_PZAP,h_PZAP=solve(tnew,N,y0z,PZAP,C1,C2,g)
    ca_PZAP=np.asarray(ca_PZAP)
    h_PZAP=np.asarray(h_PZAP)
    #plt.plot(tnew,ca,'m') #model generated
    print(f'ca0 in main',ca_PZAP[0])
    # plt.show()
    print(ca_PZAP)
    

#5. Ca Signal from the ODE model uisng the PSYK signal   ##########PSYK

    CA0s=0.247688
    y0s = [CA0, 1]
    ca_PSYK,h_PSYK=solve(tnew,N,y0s,PSYK,C1s,C2s,gs)
    ca_PSYK=np.asarray(ca_PSYK)
    h_PSYK=np.asarray(h_PSYK)
    #plt.plot(tnew,ca,'m') #model generated
    print(f'ca0 in main',ca_PSYK[0])
    # plt.show()
    print(ca_PSYK)
           

    #ca=ca_PZAP+ca_PSYK
    
    ca=ca_PSYK
    
    #plt.plot(tnew,ca_PZAP, 'r*',tnew,ca_PSYK,'b*')
    #plt.show()
    
    
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
    

    #plt.plot(tnew,ca,'c-',TT,MODEL, 'g-',TT,EXPT,'b*',time,exp_data,'r-',linewidth=2.0)

    #cmap = plt.cm.Reds
   
    
    # Set the y-axis ticks between 0.2 and 1 with a gap of 0.1
    y_ticks = np.arange(0.2, 1.0, 0.1)
    plt.yticks(y_ticks)
   # plt.plot(TT,EXPT,color=cmap(0.8),linewidth=3.0)
   # plt.plot(TT,EXPT,'skyblue',linewidth=3.0)
  
    #plt.plot(TT,EXPT,'blue',linewidth=3.0)
    #plt.plot(TT,EXPT,'blue',linewidth=3.0)
    
    
    
   # plt.plot(TT,MODEL,'k',linewidth=3.0)
   # plt.show()
    
    # Write the data to a file in the folder
    file_path = os.path.join(folder_name, "SYK_Ca.txt")
    with open(file_path, 'w') as file:
        for i in range(len(time)):
            file.write(f"{time[i]} {exp_data[i]}\n")
            
     # Write the data to a file in the folder
    file_path = os.path.join(folder_name, "SYK_Ca_fit.txt")
    with open(file_path, 'w') as file:
         for i in range(len(TT)):
             file.write(f"{TT[i]} {MODEL[i]}\n")
    
    
    

#Ca plot
    point_size = 35
    plt.scatter(time, exp_data, s=point_size, c='none', edgecolors='blue', marker='o', label='Expt')
    plt.plot(TT, MODEL, color='navy', linewidth=4, label='Model')
    #manual_x_ticks = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270]
    #manual_x_tick_labels = ['0', '30', '60', '90', '120', '150', '180', '210', '240', '270']
    manual_x_ticks = [0, 60, 120, 180, 240, 300]
    manual_x_tick_labels = ['0', '60', '120', '180', '240', '300']
    plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=oo)
    # Add minor ticks at positions 30, 90, 150, 210, 270
    minor_x_ticks = [30, 90, 150, 210, 270]
    plt.gca().set_xticks(minor_x_ticks, minor=True)
    manual_y_ticks = [0.2, 0.4, 0.6, 0.8]
    manual_y_tick_labels = ['0.2', '0.4', '0.6', '0.8']
    plt.yticks(manual_y_ticks, manual_y_tick_labels, fontsize=oo)   
    plt.tick_params(axis='both', which='major', width=2, length=8)
    
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
    plt.legend()
    
    
    print("The (SSR)**0.25 is SSR=" + str(SSR**(1/4.0)))
    
    print("K1 = "+str(K1))
    print("K2 = "+str(K2))
    print("K3 = "+str(K3))
    print("KSA = "+str(KSA))
    print("KST = "+str(KST))
    print("C1 = "+str(C1))
    print("C2 = "+str(C2))
    print("g = "+str(g))
    print("C1s = "+str(C1s))
    print("C2s = "+str(C2s))
    print("gs = "+str(gs))
    print("SSR=" + str(SSR**(1/4.0)))
    
    
    print('remove done')


    return SSR,tnew,ca,TT,MODEL,TT,EXPT