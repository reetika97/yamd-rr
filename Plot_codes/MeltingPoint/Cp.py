import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d
import os


def heat_capacity(num_atoms, file, e_plot, cp_plot):
    data = pd.read_csv(file, sep=';')

    # Extracting energy and temperature
    energy = data.iloc[:, 0]
    temperature = data.iloc[:, 1]

    e_plot.plot(temperature, energy, '-', color='darkblue',
                label=num_atoms)
    e_plot.legend(loc='upper left', fontsize='small', frameon=True, handlelength=0.5, handleheight=0.5)
    e_plot.grid(True)

    temperature = uniform_filter1d(temperature, size=7)
    energy = uniform_filter1d(energy, size=7)
    d_energy = np.diff(energy)
    d_temp = np.diff(temperature)

    # Identify the index where the slope changes significantly
    melting_point_index = np.argmin(d_temp)
    melting_point = temperature[melting_point_index + 1]  # Offset by one due to np.diff

    # Calculate latent heat (energy difference at phase change)
    latent_heat = np.where((temperature >= (melting_point - 15)) & (temperature <= (melting_point + 15)))[0]

    # Create a smooth interpolating function using UnivariateSpline
    # Differentiate the interpolating function to obtain Cp (dE/dT)
    cp1 = d_energy[:latent_heat[0]] / d_temp[:latent_heat[0]]
    t1 = temperature[:latent_heat[0]]
    cp2 = (energy[latent_heat[-1]] - energy[latent_heat[0]]) / (
            temperature[latent_heat[-1]] - temperature[latent_heat[0]])
    t2 = melting_point
    cp3 = d_energy[latent_heat[-1] + 1:] / d_temp[latent_heat[-1] + 1:]
    t3 = temperature[latent_heat[-1] + 2:]

    t1 = np.append(t1, t2)
    t1 = np.append(t1, t3)

    cp1 = np.append(cp1, cp2)
    cp1 = np.append(cp1, cp3)

    neg_index = np.where(cp1 < 0)[0]
    t1 = np.delete(t1, neg_index)
    cp1 = np.delete(cp1, neg_index)

    Cp_600 = np.mean(d_energy[np.where(np.abs(t1 - 600) < 30)[0]] / d_temp[np.where(np.abs(t1 - 600) < 30)[0]])
    Sp_cp = Cp_600 / num_atoms

    # Plot the Cp vs Temperature curve
    cp_plot.plot(t1, cp1, label=num_atoms, color='darkblue')
    cp_plot.legend(loc='upper left',  fontsize='small', frameon=True, handlelength=0.5, handleheight=0.5)
    cp_plot.grid(True)

    return [num_atoms, Cp_600, Sp_cp]


if __name__ == '__main__':

    # Setting font to Arial
    plt.rcParams["font.family"] = "Arial"

    # Create a figure
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(10, 2, figure=fig)

    # Create a list of subplots
    ax = [fig.add_subplot(gs[i, j]) for i in range(10) for j in range(2)]

    # File processing
    files = os.listdir()
    i = 0
    num_atoms = []
    cp_600 = []
    for file in files:
        if file.endswith('.csv'):
            num_atoms.append(int(file[:-8]))

    num_atoms.sort()
    for key in num_atoms:
        file = str(key) + 'EvsT.csv'
        cp_600.append(heat_capacity(key, file, ax[i], ax[i + 1]))
        i += 2

    # Plot for energy and cp.
    ax[10].set_ylabel('Total Energy (eV)', labelpad=10, fontsize=14)
    ax[11].set_ylabel('Heat Capacity $C_p$ (eV/K)', labelpad=10, fontsize=14)
    ax[18].set_xlabel('Temperature (K)', fontsize=14)
    ax[19].set_xlabel('Temperature (K)', fontsize=14)
    ax[0].set_title('Total Energy vs Temperature', fontsize=16)
    ax[1].set_title('$C_p$ vs Temperature', fontsize=16)
    plt.tight_layout(w_pad=2)
    plt.show()

    # plot Cp @ 600 and Sp @ 600:

    cp_600 = np.array(cp_600)
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(1, 2, figure=fig)

    # Create a list of subplots
    ax = [fig.add_subplot(gs[0, j]) for j in range(2)]

    ax[0].plot(cp_600[:, 0], cp_600[:, 1], '-o', color='darkblue', markerfacecolor='black')
    ax[0].set_ylabel('Heat Capacity $C_p$ (eV/K)', labelpad=10, fontsize=14)
    ax[0].set_xlabel('Number of Au atoms', labelpad=10, fontsize=14)
    ax[0].grid(True)
    ax[1].plot(cp_600[:, 0], cp_600[:, 2], '-o', color='darkblue', markerfacecolor='black')
    ax[1].set_ylabel('Specific Heat Capacity $c_p$ (eV/atom.K)', labelpad=10, fontsize=14)
    ax[1].set_xlabel('Number of Au atoms', labelpad=10, fontsize=14)
    ax[1].grid(True)
    ax[0].set_title('Heat Capacity $C_p$ @ 600K vs Cluster Size', fontsize=18)
    ax[1].set_title('Specific Heat $c_p$ @ 600K vs Cluster Size', fontsize=18)
    ax[0].text(0.02, 0.95, '(a)', transform=ax[0].transAxes, fontsize=14, va='top', ha='left')
    ax[1].text(0.02, 0.95, '(b)', transform=ax[1].transAxes, fontsize=14, va='top', ha='left')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout(w_pad=1)
    plt.show()





