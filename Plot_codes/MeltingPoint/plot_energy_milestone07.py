from matplotlib import pyplot as plt
import pandas as pd

plt.figure()

file = r'923EvsT.csv'
data = pd.read_csv(file, sep=';')
# Extracting energy and temperature
energy = data.iloc[:, 0]
temperature = data.iloc[:, 1]

plt.plot(energy[::2], temperature[::2], '-o', markersize = 5, markerfacecolor='black', color='darkblue',
         markeredgecolor='blue')

plt.xlabel('Energy (eV)')
plt.ylabel('Temperature (K)')
plt.title('Energy vs Temperature')
plt.grid(True)
plt.show()
