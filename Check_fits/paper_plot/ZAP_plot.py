#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:14:04 2024

@author: ixn004
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:44:04 2024

@author: ixn004
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data from files for ZAP only
ZAP_ca_fit = np.loadtxt("ZAP_Ca_fit.txt")
ZAP_ca = np.loadtxt("ZAP_Ca.txt")

# Extract first and second columns
TT_ZAP, MODEL_ZAP = ZAP_ca_fit [:, 0], ZAP_ca_fit [:, 1]
time_ZAP,exp_data_ZAP = ZAP_ca[:, 0], ZAP_ca[:, 1]



#Ca plot Mixed
oo=24
point_size = 35
b=4
a=2


#ZAP only
plt.scatter(time_ZAP, exp_data_ZAP, s=point_size, c='none', edgecolors='red', marker='o', label='Expt')
plt.plot(TT_ZAP, MODEL_ZAP, color='maroon', linewidth=b, label='Model')



#Ticks
manual_x_ticks = [0, 60, 120, 180, 240]
manual_x_tick_labels = ['0', '60', '120', '180', '240']
plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=oo)
# Add minor ticks at positions 30, 90, 150, 210, 270
minor_x_ticks = [30, 90, 150, 210, 270]
plt.gca().set_xticks(minor_x_ticks, minor=True)
manual_y_ticks = [0.2, 0.4, 0.6, 0.8]
manual_y_tick_labels = ['0.2','0.4','0.6','0.8']
#manual_y_ticks = [0.2, 0.3, 0.4]
#manual_y_tick_labels = ['0.2', '0.3','0.4']
plt.yticks(manual_y_ticks, manual_y_tick_labels, fontsize=oo)
plt.tick_params(axis='both', which='major', width=2, length=8)


# Create a boxplot and set the width of the boxes
# Customize the axes line width
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(a)  # Adjust the left spine's line width
# Adjust the bottom spine's line width
ax.spines['bottom'].set_linewidth(a)
# Adjust the right spine's line width
ax.spines['right'].set_linewidth(a)
ax.spines['top'].set_linewidth(a)  # Adjust the top spine's line width
#plt.legend()
plt.show()



manual_x_ticks = [0, 60, 120, 180, 240]
manual_x_tick_labels = ['0', '60', '120', '180', '240']
plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=oo)
# Add minor ticks at positions 30, 90, 150, 210, 270
minor_x_ticks = [30, 90, 150, 210, 270]
plt.gca().set_xticks(minor_x_ticks, minor=True)
manual_y_ticks = [0.1,0.2, 0.3]
manual_y_tick_labels = ['0.1','0.2','0.3']
#manual_y_ticks = [0.2, 0.3, 0.4]
#manual_y_tick_labels = ['0.2', '0.3','0.4']
plt.yticks(manual_y_ticks, manual_y_tick_labels, fontsize=oo)
plt.tick_params(axis='both', which='major', width=2, length=8)

# Create a boxplot and set the width of the boxes
# Customize the axes line width
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(a)  # Adjust the left spine's line width
# Adjust the bottom spine's line width
ax.spines['bottom'].set_linewidth(a)
# Adjust the right spine's line width
ax.spines['right'].set_linewidth(a)
ax.spines['top'].set_linewidth(a)  # Adjust the top spine's line width
#plt.legend()
plt.show()







# Load data from files for mixed
ZAP_only = np.loadtxt("ZAP_only.txt")
# Extract first and second columns
T_ZAP, ZAP_PZAP,ZAP_PSYK = ZAP_only  [:, 0], ZAP_only [:, 1],ZAP_only [:, 2]

ooo=24

#Plot ZAP  number
# manual_y_ticks = [0, 100, 200,300]
# manual_y_tick_labels = ['0','100','200','300']
#manual_y_ticks = [0.2, 0.3, 0.4]
#manual_y_tick_labels = ['0.2', '0.3','0.4']
manual_x_ticks = [0, 60, 120, 180, 240]
manual_x_tick_labels = ['0', '60', '120', '180', '240']
plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=ooo)
# Add minor ticks at positions 30, 90, 150, 210, 270
minor_x_ticks = [30, 90, 150, 210, 270]
plt.gca().set_xticks(minor_x_ticks, minor=True)

# Setting manual major ticks and labels
y_ticks = np.arange(0, 900, 150)
plt.yticks(y_ticks)

#plt.plot(T_mixed,mixed_PSYK,'black',linestyle='--', linewidth=4,label='pSYK') 
plt.plot(T_ZAP,ZAP_PZAP, color='red', linewidth=6, label='pZAP')
plt.tick_params(axis='both', which='major', labelsize=ooo)  # Adjust label size
plt.tick_params(axis='both', which='major', width=2, length=8)
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(2.0)  # Adjust the left spine's line width
ax.spines['bottom'].set_linewidth(2.0)  # Adjust the bottom spine's line width
ax.spines['right'].set_linewidth(2.0)  # Adjust the right spine's line width
ax.spines['top'].set_linewidth(2.0)  # Adjust the top spine's line width
#plt.legend()
plt.show()


# manual_y_ticks = [0, 100, 200]
# manual_y_tick_labels = ['0', '100', '200',]
# plt.yticks(manual_y_ticks, manual_y_tick_labels, fontsize=oo)
# # Add minor ticks at positions 30, 90, 150, 210, 270
# minor_y_ticks = [50, 150]
# plt.gca().set_yticks(minor_y_ticks, minor=True)

# y_ticks = np.arange(0, 300, 100)
# plt.yticks(y_ticks)




