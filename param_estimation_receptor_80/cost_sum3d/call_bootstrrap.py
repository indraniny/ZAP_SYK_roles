#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:39:54 2022

@author: ixn004


"""


from ca_bootstrap import bootstrap_ca


def main():
    m=100
    mean_t,ca_boot=bootstrap_ca(m)
    A=[]
    A=[mean_t,ca_boot]
    #print(A)
    # f = open('param_residue.dat', 'w')
    #     f.write("The optimized parameter is " + str(optimized_param) + "\n Residue = " + str(residue) +"\n Squared residue = " +str(xx)  )
    #     f.close()
    
if __name__=="__main__" :
    main()    
          