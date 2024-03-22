#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 11:32:33 2024

@author: ixn004
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data from files for mixed
mixed_ZAP_SYK = np.loadtxt("mixed_ZAP_SYK.txt")

# Extract first and second columns
T_mixed, mixed_PZAP,mixed_PSYK = mixed_ZAP_SYK [:, 0], mixed_ZAP_SYK [:, 1],mixed_ZAP_SYK [:, 2]


# Load data from files for mixed
SYK_only = np.loadtxt("SYK_only.txt")
# Extract first and second columns
T_SYK, SYK_PZAP,SYK_PSYK = SYK_only  [:, 0], SYK_only [:, 1],SYK_only [:, 2]


# Load data from files for mixed
ZAP_only = np.loadtxt("ZAP_only.txt")
# Extract first and second columns
T_ZAP, ZAP_PZAP,ZAP_PSYK = ZAP_only  [:, 0], ZAP_only [:, 1],ZAP_only [:, 2]

oo=16

#Plot ZAP SYK number
y_ticks = np.arange(0, 3000, 500)
plt.yticks(y_ticks)
plt.plot(T_mixed,mixed_PZAP, color='black', linewidth=6, label='pZAP')
plt.plot(T_mixed,mixed_PSYK,'black',linestyle='--', linewidth=4,label='pSYK') 
plt.plot(T_ZAP,ZAP_PZAP, color='red', linewidth=6, label='pZAP')
plt.tick_params(axis='both', which='major', labelsize=oo)  # Adjust label size
plt.tick_params(axis='both', which='major', width=2, length=8)
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(2.0)  # Adjust the left spine's line width
ax.spines['bottom'].set_linewidth(2.0)  # Adjust the bottom spine's line width
ax.spines['right'].set_linewidth(2.0)  # Adjust the right spine's line width
ax.spines['top'].set_linewidth(2.0)  # Adjust the top spine's line width
#plt.legend()
plt.show()


#Plot ZAP SYK number
y_ticks = np.arange(0, 7000, 2000)
plt.yticks(y_ticks)
plt.plot(T_mixed,mixed_PZAP, color='black', linewidth=6, label='pZAP')
plt.plot(T_mixed,mixed_PSYK,'black',linestyle='--', linewidth=4,label='pSYK') 
plt.plot(T_SYK,SYK_PSYK, color='blue', linestyle='--',linewidth=4, label='pSYK')
plt.tick_params(axis='both', which='major', labelsize=oo)  # Adjust label size
plt.tick_params(axis='both', which='major', width=2, length=8)
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(2.0)  # Adjust the left spine's line width
ax.spines['bottom'].set_linewidth(2.0)  # Adjust the bottom spine's line width
ax.spines['right'].set_linewidth(2.0)  # Adjust the right spine's line width
ax.spines['top'].set_linewidth(2.0)  # Adjust the top spine's line width
#plt.legend()
plt.show()
