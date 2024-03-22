import os
import numpy as np
import matplotlib.pyplot as plt

# Define the folder names

#folder_names = ["ZAP_HPC_output_0", "ZAP_HPC_output_1", "ZAP_HPC_output_2", "ZAP_HPC_output_3", "ZAP_HPC_output_4"]
#folder_names = ["SYK_HPC_output_0", "SYK_HPC_output_1", "SYK_HPC_output_2", "SYK_HPC_output_3", "SYK_HPC_output_4"]
folder_names = ["WT_HPC_output_0", "WT_HPC_output_1", "WT_HPC_output_2", "WT_HPC_output_3", "WT_HPC_output_4"]

# Initialize an empty list to store data from each folder
data_list = []

# Iterate over each folder
for folder_name in folder_names:
    #file_path = os.path.join(folder_name, "Model_SYK.gdat")
    file_path = os.path.join(folder_name, "Model_WT.gdat")
    #file_path = os.path.join(folder_name, "Model_ZAP.gdat")
    # Check if the file exists
    if os.path.exists(file_path):
        # Read the data from the file
        with open(file_path, 'r') as file:
            # Assuming the data is space-separated
            data = np.loadtxt(file)
            data_list.append(data)
    else:
        print(f"The file 'Model_WT.gdat' does not exist in the folder {folder_name}.")

# Check if any data was read
if data_list:
    # Take the average across all arrays in the list
    average_data = np.mean(data_list, axis=0)
    
    # Save the average data into a text file
    #output_file_path = "average_data_ZAP.txt"#************************************************************
    output_file_path = "average_data_WT.txt"#************************************************************
    #output_file_path = "average_data_SYK.txt"#************************************************************
    
    np.savetxt(output_file_path, average_data, fmt='%f', delimiter='\t')
    print(f"Average data saved to '{output_file_path}' successfully.")
else:
    print("No data was read from any folder.")
    

# Molecules pSYK_total A(State~pSYK!?)
# Molecules pSYK_free A(State~pSYK)
# Molecules pSYK_bound A(State~pSYK!+)

# Molecules PZAP_total A(State~PZAP!?)
# Molecules PZAP_free A(State~PZAP)
# Molecules PZAP_bound A(State~PZAP!+)

# Molecules tot_bound_ZAP A(State~UZAP!+) A(State~PZAP!+)
# Molecules tot_bound_SYK A(State~uSYK!+) A(State~pSYK!+)

# Molecules bound_phosph_zeta Zeta(receptor!+,ITAM1~PP!?) Zeta(receptor!+,ITAM2~PP!?) Zeta(receptor!+,ITAM3~PP!?) Zeta(receptor!+,ITAM4~PP!?) Zeta(receptor!+,ITAM5~PP!?) Zeta(receptor!+,ITAM6~PP!?)
# Molecules open_U_zeta Zeta(receptor!+,ITAM1~U) Zeta(receptor!+,ITAM2~U) Zeta(receptor!+,ITAM3~U) Zeta(receptor!+,ITAM4~U) Zeta(receptor!+,ITAM5~U) Zeta(receptor!+,ITAM6~U)
# Molecules bound_U_zeta Zeta(receptor!+,ITAM1~U!?) Zeta(receptor!+,ITAM2~U!?) Zeta(receptor!+,ITAM3~U!?) Zeta(receptor!+,ITAM4~U!?) Zeta(receptor!+,ITAM5~U!?) Zeta(receptor!+,ITAM6~U!?)
# Molecules total_phosph_zeta Zeta(ITAM1~PP!?) Zeta(ITAM2~PP!?) Zeta(ITAM3~PP!?) Zeta(ITAM4~PP!?) Zeta(ITAM5~PP!?) Zeta(ITAM6~PP!?)

# Molecules tot_bound_SHP A(State~SHP!+) 
# Molecules lig_receptor_complex CD16(lig!1).Ligand(Site!1)
# Molecules Receptor CD16(lig!+) CD16(lig)
# Molecules zeta Zeta(receptor!+) Zeta(receptor)
# Molecules Lig Ligand(Site!+) Ligand(Site)


time = average_data[:, 0]

total_psyk= average_data[:, 1]
free_psyk=average_data[:, 2]
bound_psyk=average_data[:, 3]

total_pzap= average_data[:, 4]
free_pzap=average_data[:, 5]
bound_pzap=average_data[:, 6]

bound_zap=average_data[:, 7]
bound_syk= average_data[:, 8]


bound_pITAM=average_data[:, 9]
open_uITAM=average_data[:, 10]
bound_uITAM=average_data[:, 11]
total_pITAM=average_data[:, 12]


bound_shp=average_data[:, 13]
lig_recep_complex=average_data[:, 14]
recep=average_data[:, 15]
zeta=average_data[:, 16]
lig=average_data[:, 17]


n=5
oo=30
# Plot each column with time
# plt.plot(time, total_psyk, color='red', label='Total Psyk', linewidth=n)
# plt.plot(time, free_psyk, color='orange', label='Free Psyk', linewidth=n)
# plt.plot(time, bound_psyk, color='cyan', label='Bound Psyk', linewidth=n)
# plt.plot(time, total_pzap, color='green', label='Total Pzap', linewidth=n)
# plt.plot(time, free_pzap, color='blue', label='Free Pzap', linewidth=n)
# plt.plot(time, bound_pzap, color='indigo', label='Bound Pzap', linewidth=n)

#print(bound_shp)
QQ=bound_zap+bound_syk

plt.plot(time, bound_zap, color='violet',label='Bound Zap', linewidth=n)
plt.plot(time, bound_syk, color='pink', label='Bound Syk', linewidth=n)
plt.plot(time, bound_shp, color='grey', label='bound SHP', linewidth=n-3)
plt.plot(time, bound_pITAM, color='magenta', label='bound pITAM', linewidth=n)
#plt.plot(time, QQ, color='black', label='bound ZAP+SYK', linewidth=n-3)

# y_ticks = np.arange(0, 5000, 1500)
# plt.yticks(y_ticks)
# minor_y_ticks = [500, 1500]
# plt.gca().set_yticks(minor_y_ticks, minor=True)


#plt.plot(time, total_pITAM, color='pink', label='total pITAM', linewidth=n)
y_ticks = np.arange(0, 7000, 2000)
plt.yticks(y_ticks)
# minor_y_ticks = [500, 1500]
# plt.gca().set_yticks(minor_y_ticks, minor=True)


#Ticks
manual_x_ticks = [0, 60, 120, 180, 240]
manual_x_tick_labels = ['0', '60', '120', '180', '240']
plt.xticks(manual_x_ticks, manual_x_tick_labels, fontsize=oo)
# Add minor ticks at positions 30, 90, 150, 210, 270
minor_x_ticks = [30, 90, 150, 210, 270]
plt.gca().set_xticks(minor_x_ticks, minor=True)
#y_ticks = np.arange(0, 1200, 300)
# manual_y_ticks = [0,1600,3200]
# manual_y_tick_labels = ['0', '1600','3200']
# plt.yticks(manual_y_ticks, manual_y_tick_labels, fontsize=oo)
# #Add minor ticks at positions 30, 90, 150, 210, 270
# minor_y_ticks = [500,1500,2500]
# plt.gca().set_yticks(minor_y_ticks, minor=True)
plt.tick_params(which='major', width=4)

plt.tick_params(which='minor', width=4)

# Add labels and legend
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(2.0)  # Adjust the left spine's line width
# Adjust the bottom spine's line width
ax.spines['bottom'].set_linewidth(2.0)
# Adjust the right spine's line width
ax.spines['right'].set_linewidth(2.0)
ax.spines['top'].set_linewidth(2.0)  # Adjust the top spine's line width
# plt.xlabel('Time')
# plt.ylabel('Numbers')
plt.legend()
    
# Show plot
plt.show()

