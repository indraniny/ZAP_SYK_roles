#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 20:23:03 2024

@author: ixn004
"""

import numpy as np
import matplotlib.pyplot as plt

data_WT = np.loadtxt("average_data_WT.txt")
data_ZAP= np.loadtxt("average_data_ZAP.txt")
data_SYK= np.loadtxt('average_data_SYK.txt')


# total_psyk= average_data[:, 1]
# free_psyk=average_data[:, 2]
# bound_psyk=average_data[:, 3]

# total_pzap= average_data[:, 4]
# free_pzap=average_data[:, 5]
# bound_pzap=average_data[:, 6]

# bound_zap=average_data[:, 7]
# bound_syk= average_data[:, 8]


# bound_pITAM=average_data[:, 9]
# open_uITAM=average_data[:, 10]
# bound_uITAM=average_data[:, 11]
# total_pITAM=average_data[:, 12]


# bound_shp=average_data[:, 13]
# lig_recep_complex=average_data[:, 14]
# recep=average_data[:, 15]
# zeta=average_data[:, 16]
# lig=average_data[:, 17]


time_WT = data_WT[:, 0]
bound_zap_WT=data_WT[:, 4]
bound_syk_WT= data_WT[:, 1]


time_ZAP = data_ZAP[:, 0]
bound_zap_ZAP=data_ZAP[:, 4]
bound_syk_ZAP= data_ZAP[:, 1]

time_SYK = data_SYK[:, 0]
bound_zap_SYK=data_SYK[:, 4]
bound_syk_SYK= data_SYK[:, 1]




# time_WT = data_WT[:, 0]
# bound_zap_WT=data_WT[:, 7]
# bound_syk_WT= data_WT[:, 8]


# time_ZAP = data_ZAP[:, 0]
# bound_zap_ZAP=data_ZAP[:, 7]
# bound_syk_ZAP= data_ZAP[:, 8]

# time_SYK = data_SYK[:, 0]
# bound_zap_SYK=data_SYK[:, 7]
# bound_syk_SYK= data_SYK[:, 8]


#Ticks
oo=30
point_size = 35
b=3
a=3 #boxwidth

#ZAP
plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')


manual_x_ticks = [0, 60, 120, 180, 240]
manual_x_tick_labels = ['0', '60', '120', '180', '240']
plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=oo)
# Add minor ticks at positions 30, 90, 150, 210, 270
minor_x_ticks = [30, 90, 150, 210, 270]
plt.gca().set_xticks(minor_x_ticks, minor=True)

# Setting manual major ticks and labels
y_ticks = np.arange(0, 1000, 150)
plt.yticks(y_ticks)

# #WT
# plt.plot(time_WT,bound_zap_WT, color='violet', linewidth=6, label='bound ZAP')
# plt.plot(time_WT,bound_syk_WT, color='pink', linewidth=6, label='bound SYK')


#ZAP_only
# plt.plot(time_ZAP,bound_zap_ZAP, color='violet', linewidth=6, label='bound ZAP')
# plt.plot(time_ZAP,bound_syk_ZAP, color='pink', linewidth=6, label='bound SYK')

# #SYK_only
# plt.plot(time_SYK,bound_zap_SYK, color='violet', linewidth=6, label='bound ZAP')
# plt.plot(time_SYK,bound_syk_SYK, color='pink', linewidth=6, label='bound SYK')



#WT vs SYK
# plt.plot(time_WT,bound_syk_WT, color='slategrey', linewidth=6,linestyle='--', label='bound SYK_WT')
# plt.plot(time_SYK,bound_syk_SYK, color='royalblue', linewidth=6,linestyle='--', label='bound SYK_SYK')

# plt.plot(time_WT,bound_syk_WT, color='black', linewidth=6,linestyle='--', label='bound SYK_WT')
# plt.plot(time_SYK,bound_syk_SYK, color='blue', linewidth=6,linestyle='--', label='bound SYK_SYK')

#WT vs ZAP
# plt.plot(time_WT,bound_zap_WT, color='black', linewidth=6, label='bound ZAP_WT')
# plt.plot(time_ZAP,bound_zap_ZAP, color='red', linewidth=6, label='bound ZAP_WT')

plt.plot(time_WT,bound_zap_WT, color='slategrey', linewidth=6, label='bound ZAP_WT')
plt.plot(time_ZAP,bound_zap_ZAP, color='orangered', linewidth=6, label='bound ZAP_WT')


plt.tick_params(axis='both', which='major', labelsize=oo)  # Adjust label size
plt.tick_params(axis='both', which='major', width=2, length=8)
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(a)  # Adjust the left spine's line width
ax.spines['bottom'].set_linewidth(a)  # Adjust the bottom spine's line width
ax.spines['right'].set_linewidth(a)  # Adjust the right spine's line width
ax.spines['top'].set_linewidth(a)  # Adjust the top spine's line width
#plt.legend()
plt.show()

