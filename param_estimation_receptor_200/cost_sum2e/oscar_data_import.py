#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 23:16:34 2022

@author: ixn004
"""

import pandas as pd
import numpy as np
import scipy as sp


def expt_import(n):
     #start by importing your data
    data=pd.read_csv('oscar_ca.dat',sep="\t", comment='#', header=None)
    #if the data has a time column:
    data=np.asarray(data) 
    #print(type(data))
    time=data[:,0]
    #print(len(time))
    #print(time)
    exp_data=data[:,1]
    #plt.plot(time,exp_data,'g--')
    return time,exp_data