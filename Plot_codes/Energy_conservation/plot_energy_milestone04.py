import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

plt.figure()

file = r'0.001_Etot_file.csv'
Energy = pd.read_csv(file, sep=';')
t = np.array(range(len(Energy)))

# Setting font to Arial
plt.rcParams["font.family"] = "Arial"

# Plotting the lines with different formats
plt.plot(t, Energy['Epot'], label='Epot', linewidth=2, linestyle='--', color='blue')
plt.plot(t, Energy['Ekin'], label='Ekin', linewidth=2, linestyle='-.', color='purple')
plt.plot(t, Energy['Etot'], label='Etot', linewidth=2, linestyle='-', color='red')

# Adding text labels directly on the lines
plt.text(t[80], Energy['Epot'][80]+15, 'Epot', fontsize=12, color='blue')
plt.text(t[80], Energy['Ekin'][80]+15, 'Ekin', fontsize=12, color='purple')
plt.text(t[80], Energy['Etot'][80]+15, 'Etot', fontsize=12, color='red')

plt.grid(True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Time (LJ units)', fontsize=14)
plt.ylabel('Energy (LJ units)', fontsize=14)
plt.title('Energy Conservation for MD', fontsize=18)
plt.show()
