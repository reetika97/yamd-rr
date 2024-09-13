import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

colors = ['orange', 'darkblue', 'darkgreen', 'purple', 'brown']
linestyles = ['-', '--', ':', '-.', ':']

if __name__ == '__main__':
    files = os.listdir('strain_stress_files')
    Temp = []
    plt.rcParams["font.family"] = "Arial"
    fig = plt.figure()

    for file in files:
        if file.startswith('3050_1.0'):
            Temp.append(float(file.split('_')[2]))

    Temp.sort()
    idx = 0
    for key in Temp:
        if key in [0.0, 5.0, 10.0, 20.0, 50.0]:
            file = f'3050_1.0_{key:.1f}_svf.csv'
            print(file)
            file_path = 'strain_stress_files/' + file
            stress = pd.read_csv(file_path, sep=';')

            strain = stress['strain']
            # strain = uniform_filter1d(strain, size=3)
            stress = (stress['force_lg'] / stress['area']) * 160.218  # Convert to GPa
            # stress = uniform_filter1d(stress, size=3)

            # Use unique color and linestyle for each plot
            plt.plot(
                strain, stress,
                marker='o', markersize='2',
                label=f'{key:.1f} K',
                color=colors[idx % len(colors)],
                linestyle=linestyles[idx % len(linestyles)],
            )
            # plt.text(strain.iloc[-1]-0.002, stress.iloc[-1]+0.02, f'{key} atoms',  color=colors[idx % len(colors)])
            idx += 1

    plt.grid(True)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Strain $\epsilon$', fontsize=16, labelpad=10)
    plt.ylabel('Stress $\sigma$ (GPa)', fontsize=16, labelpad=10)
    plt.title('Stress vs Strain for varying Temperatures', fontsize=20)
    plt.legend(loc='upper right', fontsize=12, frameon=True)
    plt.show()
