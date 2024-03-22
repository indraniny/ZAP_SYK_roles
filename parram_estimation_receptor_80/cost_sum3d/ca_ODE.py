#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 22:17:27 2022

@author: ixn004
"""

import numpy as np
from scipy.integrate import odeint

def calcium(y, t, PZAP,C1,C2,g):
    #be=0.02 # leak term from cytosol
    # be=0
    # b=0.111 # maximum Ca flow from Cyto to ER
    # k1=0.7 #
    # k2=0.7 #
    # s=(k2*k2) #
    # #mu=0.3  #bifurcation parameter- depends on the IP3 concentration, which is assumed to be constant
    # f=0.002 # 
    
    # #(th=1500,g=0.006,kf=6.6,v=0.05)
  
    # kf=6.66 #k effective
    # th=1500.0 #
    # kg=10 #k_gamma
    # g=0.007
    # #g=0.005 #
    # v=0.05
    
    
    
    #z=602 #from uM unit to molecules/mu^3 unit
   
    be=0
    b=0.111 # maximum Ca flow from Cyto to ER  
    k1=0.7 # In Atri 1993, has unit of uM 
    k2=0.7 # In Atri 1993, has unit of uM 
    s=(k2*k2) #
    
   
   
    
    # #mu=0.3  #bifurcation parameter- depends on the IP3 concentration, which is assumed to be constant
    # f=0.002 # 
    
    # #(th=1500,g=0.006,kf=6.6,v=0.05)
      
    
    #th=11.0 #
    #kg=10 #k_gamma
    
    ##g=0.005 #
    
    #g=0.06
    #v=0.00008
    #wt
    #kf=wt*z #k effective
      
    
    ca, h = y
    #pZAP is in the uM unit
    
    #dcadt=(kf*h*v*mu*PZAP*(((b*k1)+ca)/(k1+ca)))-((g*ca)/(kg+ca))+(be)
    dcadt=(C1*h*PZAP*(((b*k1)+ca)/(k1+ca)))-(g*ca)+(be)
    dhdt=(C2*PZAP)*((s/(s+ca ** 2))-h)
    
    # np.asanyarray(Y_smooth,dtype=np.float64)
    # print(Y_smooth)
    # print(len(Y_smooth))
    # type(Y_smooth)
    


    #dcadt=(kf*mu*h*(1/(c+np.exp(-b1*(t ** m))))*((b*k1+ca)/(k1+ca)))-(g*f*ca)+(be)
    #dhdt=(1.0/(th*(c+np.exp(-b1*(t ** m)))))*((s/(s+ca ** 2))-h)

    # print(dcadt)
    #print(type(dcadt))
    #print(len(dcadt))
    # dydt(1)=(kf*mu*y(2)*(1/(c+exp(-b1*(t.^m))))*((b*k1+y(1))/(k1+y(1))))-g*f*y(1)+be;
    # dydt(2)=(1.0/(th*(c+exp(-b1*(t.^m)))))*((s/(s+(y(1)^2)))-y(2));
    dydt = [dcadt,dhdt]
    return dydt

def solve(tnew,N,y0,PZAP,C1,C2,g):
# store solution
    #empty_like_sets_array_of_same_shape_and_size
    ca = np.empty_like(tnew)
    h = np.empty_like(tnew)
    # record initial conditions
    ca[0] = y0[0]
    h[0] = y0[1]
    #PZAP=np.ones(len(PZAP))

    # solve ODE
    for i in range(1,N):
        # span for next time step
        tspan = [tnew[i-1],tnew[i]]
        # solve for next step
        sol = odeint(calcium,y0,tspan,args=(PZAP[i],C1,C2,g))
        # store solution for plotting
        ca[i] = sol[1][0]
        h[i] = sol[1][1]
        # next initial condition
        y0 = sol[1]
        
    return ca,h    