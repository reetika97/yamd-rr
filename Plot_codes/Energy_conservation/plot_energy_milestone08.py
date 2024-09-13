import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Reading the data from the CSV file
file = r'Energy_MPI.csv'
Energy = pd.read_csv(file, sep=';')
t = np.array(range(len(Energy))) * 100

# Setting font to Arial
plt.rcParams["font.family"] = "Arial"

# Create subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plotting Ekin on the first subplot
axs[0].plot(t, Energy['Ekin'], label='Ekin', linewidth=2, linestyle='-.', color='purple')
axs[0].set_ylabel('Energy (eV)', fontsize=14, labelpad=10)  # Increased labelpad
axs[0].set_title('Energy Conservation for MD', fontsize=18)

# Adding text labels directly on the lines for Ekin
axs[0].text(t[80], Energy['Ekin'][80]-1.7, 'Ekin', fontsize=12, color='purple')

# Plotting Epot and Etot on the second subplot
axs[1].plot(t, Energy['Epot'], label='Epot', linewidth=2, linestyle='--', color='blue')
axs[1].plot(t, Energy['Etot'], label='Etot', linewidth=2, linestyle='-', color='red')
axs[1].set_xlabel('Time (fs)', fontsize=14)
axs[1].set_ylabel('Energy (eV)', fontsize=14, labelpad=10)  # Increased labelpad


# Adding text labels directly on the lines for Epot and Etot
axs[1].text(t[80], Energy['Epot'][80]+1.8, 'Epot', fontsize=12, color='blue')
axs[1].text(t[80], Energy['Etot'][80]-1.2, 'Etot', fontsize=12, color='red')

axs[0].grid(True)
axs[1].grid(True)
# Align the y-axis labels
fig.align_ylabels(axs)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Adjust layout for better spacing
plt.tight_layout()

# Display the plot
plt.show()