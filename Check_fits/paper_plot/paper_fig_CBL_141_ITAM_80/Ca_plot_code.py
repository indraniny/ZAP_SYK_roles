#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:44:04 2024

@author: ixn004
"""
from sklearn.metrics import r2_score 
import numpy as np
import matplotlib.pyplot as plt

def r_squared(X, Y):
    mean_X = np.mean(X)
    total_sum_squares = np.sum((X - mean_X)**2) #expt
    residual_sum_squares = np.sum((Y - X)**2)
    r2 = 1 - (residual_sum_squares / total_sum_squares)
    return r2

# Load data from files for mixed
mixed_ca_pred = np.loadtxt("mixed_Ca_pred.txt")
mixed_ca = np.loadtxt("mixed_Ca.txt")
# Load data from the text file
data_ca= np.loadtxt('ca_comp.txt')

# Extract first and second columns
TT, MODEL = mixed_ca_pred[:, 0], mixed_ca_pred[:, 1]
time,exp_data = mixed_ca[:, 0], mixed_ca[:, 1]


# Load data from files for SYK only
SYK_ca_fit = np.loadtxt("SYK_Ca_fit.txt")
SYK_ca = np.loadtxt("SYK_Ca.txt")

# Extract first and second columns
TT_SYK, MODEL_SYK = SYK_ca_fit [:, 0], SYK_ca_fit [:, 1]
time_SYK,exp_data_SYK = SYK_ca[:, 0], SYK_ca[:, 1]


# Load data from files for ZAP only
ZAP_ca_fit = np.loadtxt("ZAP_Ca_fit.txt")
ZAP_ca = np.loadtxt("ZAP_Ca.txt")

# Extract first and second columns
TT_ZAP, MODEL_ZAP = ZAP_ca_fit [:, 0], ZAP_ca_fit [:, 1]
time_ZAP,exp_data_ZAP = ZAP_ca[:, 0], ZAP_ca[:, 1]


##Ca plot for WT
TTT, ca_PZAP, ca_PSYK = data_ca [:, 0], data_ca[:, 1], data_ca[:, 2]



#Ca plot Mixed
oo=30
point_size = 35
b=3
a=3 #boxwidth
# # #mixed
# plt.scatter(time, exp_data, s=point_size, c='none',edgecolors='black', marker='o', label='Expt')
# plt.plot(TT, MODEL, color='black', linewidth=b, label='Model')


interpolated_model_data = np.interp(time, TT, MODEL)
# plt.scatter(time,exp_data, s=point_size, c='none',edgecolors='black', label='Expt')
# plt.plot(time, interpolated_model_data, color='black', linewidth=b, label='Model')

# # Calculate R-squared
r2 = r_squared(exp_data, interpolated_model_data)
print("R-squared:", r2)
r2 = r2_score(exp_data, interpolated_model_data)
print("R-squared_python:", r2)
print(len(exp_data))

# #SYK only
# plt.scatter(time_SYK, exp_data_SYK, s=point_size, c='none', edgecolors='cornflowerblue', marker='o', label='Expt')
# #plt.scatter(time_SYK, exp_data_SYK, s=point_size, c='none', edgecolors='blue', marker='o', label='Expt')
# plt.plot(TT_SYK, MODEL_SYK, color='blue', linewidth=b, label='Model')

# # #ZAP only
plt.scatter(time_ZAP, exp_data_ZAP, s=point_size, c='none', edgecolors='red', marker='o', label='Expt')
plt.plot(TT_ZAP, MODEL_ZAP, color='maroon', linewidth=b, label='Model')



#Ticks
plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

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


#Plot R^2
print(exp_data)
print(interpolated_model_data)

#plt.plot(exp_data,interpolated_model_data, 'ks')
# plt.scatter(time,exp_data, s=point_size, c='none',edgecolors='black', label='Expt')
# plt.plot(time, interpolated_model_data, color='black', linewidth=b, label='Model')
# Calculate R-squared
r2 = r_squared(exp_data, interpolated_model_data)
print("R-squared:", r2)


plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

# x_ticks = np.arange(0, 0.5, 0.15)
# plt.yticks(x_ticks)



manual_y_ticks = [ 0.3, 0.4,0.5]
manual_y_tick_labels = ['0.3', '0.4','0.5']
plt.yticks(manual_y_ticks, manual_y_tick_labels, fontsize=oo)

# Customize the axes line width
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(a)  # Adjust the left spine's line width
# Adjust the bottom spine's line width
ax.spines['bottom'].set_linewidth(a)
# Adjust the right spine's line width
ax.spines['right'].set_linewidth(a)
ax.spines['top'].set_linewidth(a)  # Adjust the top spine's line width
show()


#Ca only

b=4
a=3 #boxwidth

plt.plot(TTT, ca_PZAP, color='maroon', linewidth=6, label='Model')
plt.plot(TTT, ca_PSYK, color='blue', linewidth=6, label='Model')
#Ticks
plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

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

ooo=24

#Plot ZAP  number
# manual_y_ticks = [0, 100, 200,300]
# manual_y_tick_labels = ['0','100','200','300']
#manual_y_ticks = [0.2, 0.3, 0.4]
#manual_y_tick_labels = ['0.2', '0.3','0.4']
#Ticks
plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')


manual_x_ticks = [0, 60, 120, 180, 240]
manual_x_tick_labels = ['0', '60', '120', '180', '240']
plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=ooo)
# Add minor ticks at positions 30, 90, 150, 210, 270
minor_x_ticks = [30, 90, 150, 210, 270]
plt.gca().set_xticks(minor_x_ticks, minor=True)

# Setting manual major ticks and labels
y_ticks = np.arange(0, 600, 150)
plt.yticks(y_ticks)

plt.plot(T_mixed,mixed_PZAP, color='black', linewidth=6, label='pZAP')
#plt.plot(T_mixed,mixed_PSYK,'black',linestyle='--', linewidth=4,label='pSYK') 
plt.plot(T_ZAP,ZAP_PZAP, color='red', linewidth=6, label='pZAP')
plt.tick_params(axis='both', which='major', labelsize=ooo)  # Adjust label size
plt.tick_params(axis='both', which='major', width=2, length=8)
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(a)  # Adjust the left spine's line width
ax.spines['bottom'].set_linewidth(a)  # Adjust the bottom spine's line width
ax.spines['right'].set_linewidth(a)  # Adjust the right spine's line width
ax.spines['top'].set_linewidth(a)  # Adjust the top spine's line width
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





#Plot SYK number
plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

manual_x_ticks = [0, 60, 120, 180, 240]
manual_x_tick_labels = ['0', '60', '120', '180', '240']
plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=ooo)
# Add minor ticks at positions 30, 90, 150, 210, 270
minor_x_ticks = [30, 90, 150, 210, 270]
plt.gca().set_xticks(minor_x_ticks, minor=True)


y_ticks = np.arange(0, 1600, 500)
plt.yticks(y_ticks)
#plt.plot(T_mixed,mixed_PZAP, color='black', linewidth=6, label='pZAP')
plt.plot(T_mixed,mixed_PSYK,'black',linestyle='--', linewidth=6,label='pSYK') 
plt.plot(T_SYK,SYK_PSYK, color='blue', linestyle='--',linewidth=6, label='pSYK')
plt.tick_params(axis='both', which='major', labelsize=ooo)  # Adjust label size
plt.tick_params(axis='both', which='major', width=2, length=8)
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(a)  # Adjust the left spine's line width
ax.spines['bottom'].set_linewidth(a)  # Adjust the bottom spine's line width
ax.spines['right'].set_linewidth(a)  # Adjust the right spine's line width
ax.spines['top'].set_linewidth(a)  # Adjust the top spine's line width
#plt.legend()
plt.show()



# Plotting from loaded data
plt.plot(data[:, 0], data[:, 1], label='ca_PZAP')
plt.plot(data[:, 0], data[:, 2], label='ca_PSYK')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()
plt.show()