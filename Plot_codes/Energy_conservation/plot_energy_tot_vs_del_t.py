import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

# Set figure size
plt.figure(figsize=(10, 6))

# File paths
file1 = r'0.0001_Etot_file.csv'
file2 = r'0.001_Etot_file.csv'
file3 = r'0.01_Etot_file.csv'
file4 = r'0.02_Etot_file.csv'
file5 = r'0.002_Etot_file.csv'
file6 = r'0.005_Etot_file.csv'

# Load the data
Energy_0_0001 = pd.read_csv(file1, sep=';')
Energy_0_001 = pd.read_csv(file2, sep=';')
Energy_0_01 = pd.read_csv(file3, sep=';')
Energy_0_02 = pd.read_csv(file4, sep=';')
Energy_0_002 = pd.read_csv(file5, sep=';')
Energy_0_005 = pd.read_csv(file6, sep=';')

# Time steps (assuming they are all of the same length)
t = np.array(range(len(Energy_0_001)))

# Setting font to Arial
plt.rcParams["font.family"] = "Arial"

# Extract 'Etot' columns
tot_0_0001 = Energy_0_0001['Etot']
tot_0_001 = Energy_0_001['Etot']
tot_0_01 = Energy_0_01['Etot']
tot_0_02 = Energy_0_02['Etot']
tot_0_002 = Energy_0_002['Etot']
tot_0_005 = Energy_0_005['Etot']

# Plotting the lines with different formats and reordering in decreasing Δt
plt.plot(t, tot_0_02, label='$\Delta t$ = 0.02', linestyle=':', color='purple')
plt.plot(t, tot_0_01, label='$\Delta t$ = 0.01', linestyle='-.', color='red')
plt.plot(t, tot_0_005, label='$\Delta t$ = 0.005', linestyle='--', linewidth=2, color='brown')
plt.plot(t, tot_0_002, label='$\Delta t$ = 0.002', linestyle='-', color='orange')
plt.plot(t, tot_0_001, label='$\Delta t$ = 0.001', linestyle='--', color='green')
plt.plot(t, tot_0_0001, label='$\Delta t$ = 0.0001', linestyle='-', linewidth=2, color='blue')
plt.grid(True)
# Add titles and labels
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Variation of Total Energy over time', fontname='Arial', fontsize=18)
plt.xlabel('Time (LJ-units)', fontname='Arial', fontsize=14)
plt.ylabel('Total Energy (LJ-units)', fontname='Arial', fontsize=14)

# Adjust the legend position to the lower right side
plt.legend(loc='lower right', bbox_to_anchor=(1, 0.35))

# Show the plot
plt.tight_layout()
plt.show()

