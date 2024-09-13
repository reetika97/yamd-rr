import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.ndimage import uniform_filter1d


def melting_point(file):
    plt.rcParams["font.family"] = "Arial"
    num_atoms = int(file[:-8])
    # Data (Etot, T)
    data = pd.read_csv(file, sep=';')

    # Extracting energy and temperature
    energy = data.iloc[::2, 0]
    temperature = data.iloc[::2, 1]

    # Plot energy vs temperature
    plt.plot(energy, temperature, '-o', markersize=5, markerfacecolor='black', color='darkblue',
             markeredgecolor='blue')
    plt.ylabel('Temperature (K)', fontsize=14)
    plt.xlabel('Total Energy (eV)', fontsize=14)
    plt.title('Total Energy vs Temperature', fontsize=18)
    plt.grid(True)

    # Finding the melting point (where the slope changes significantly)
    # This can be done by looking for a plateau or a sudden change in the slope.

    # Calculate derivative to find slope change
    temperature = uniform_filter1d(temperature, size=3)
    energy = uniform_filter1d(energy, size=3)
    d_energy = np.diff(energy)
    d_temp = np.diff(temperature)

    # Identify the index where the slope changes significantly
    melting_point_index = np.argmin(d_temp)
    melting_point = temperature[melting_point_index + 1]  # Offset by one due to np.diff

    # Calculate latent heat (energy difference at phase change)
    latent_heat_index = np.where((temperature >= (melting_point - 15)) & (temperature <= (melting_point + 15)))[0]

    plt.axvline(energy[latent_heat_index[0]], color='darkblue')
    plt.axvline(energy[latent_heat_index[-1]], color='darkblue')
    plt.axhline(melting_point, color='darkblue')

    plt.text(energy[0]-2*(d_energy[0]), melting_point + 4, f"Melting Point = {melting_point:.2f} K", fontsize=12)
    plt.plot([energy[latent_heat_index[0]], energy[latent_heat_index[-1]]], [temperature[-4], temperature[-4]], ':',
             color='darkblue')

    latent_heat = abs(energy[latent_heat_index[0]] - energy[latent_heat_index[-1]])
    plt.text(energy[latent_heat_index[0]]+0.25*(d_energy[0]), temperature[-4]+4, f"Latent Heat = {latent_heat:.2f} eV",
             fontsize=12)

    print(num_atoms)
    print(f"Estimated Melting Point: {melting_point} K")
    print(f"Estimated Latent Heat: {latent_heat} eV")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

    return [num_atoms, melting_point, latent_heat]


if __name__ == '__main__':
    files = os.listdir()
    mp = []
    num_atoms = []
    for file in files:
        if file.endswith('.csv'):
            num_atoms.append(int(file[:-8]))

    num_atoms.sort()
    for key in num_atoms:
        file = str(key) + 'EvsT.csv'
        mp.append(melting_point(file))

    plt.rcParams["font.family"] = "Arial"
    mp = np.array(mp)
    fig = plt.figure(figsize=(12, 10))
    plt.plot(mp[:, 0], mp[:, 1], '-o', color='darkblue', markerfacecolor='black')
    plt.ylabel('Melting Point (K)', labelpad=10)
    plt.xlabel('Number of Au atoms', labelpad=10)
    plt.title('Melting Point of Au vs Cluster Size')
    plt.grid(True)
    plt.show()

    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(1, 2, figure=fig)

    # Create a list of subplots
    ax = [fig.add_subplot(gs[0, j]) for j in range(2)]

    ax[0].plot(mp[:, 0], mp[:, 2], '-o', color='darkblue', markerfacecolor='black')
    ax[0].set_ylabel('Latent Heat (eV)', labelpad=10, fontsize=14)
    ax[0].set_xlabel('Number of Au atoms', labelpad=10, fontsize=14)
    ax[0].grid(True)
    ax[1].plot(mp[:, 0], mp[:, 2] / mp[:, 0], '-o', color='darkblue', markerfacecolor='black')
    ax[1].set_ylabel('Latent Heat (eV/atom)', labelpad=10, fontsize=14)
    ax[1].set_xlabel('Number of Au atoms', labelpad=10, fontsize=14)
    ax[1].grid(True)
    ax[0].set_title('Latent Heat vs Cluster Size', fontsize=18)
    ax[1].set_title('Latent Heat per atom vs Cluster Size', fontsize=16)
    ax[0].text(0.02, 0.95, '(a)', transform=ax[0].transAxes, fontsize=14, va='top', ha='left')
    ax[1].text(0.02, 0.95, '(b)', transform=ax[1].transAxes, fontsize=14, va='top', ha='left')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout(w_pad=1)
    plt.show()
