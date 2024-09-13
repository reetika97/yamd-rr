import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d


colors = ['orange', 'darkblue', 'darkgreen', 'purple', 'brown']
linestyles = ['-', '--', ':', '-.', ':']

if __name__ == '__main__':
    files = os.listdir('strain_stress_files')
    num_atoms = []
    plt.rcParams["font.family"] = "Arial"
    fig = plt.figure()


    for file in files:
        if file.endswith('1.0_0.0_svf.csv'):
            num_atoms.append(int(file.split('_')[0]))

    num_atoms.sort()
    idx=0
    for key in num_atoms:
        if key in [3050, 51500, 7680, 9540, 10680]:
            file = str(key) + '_1.0_0.0_svf.csv'
            file_path = 'strain_stress_files/' + file
            stress = pd.read_csv(file_path, sep=';')

            strain = stress['strain']
            #strain = uniform_filter1d(strain, size=3)
            stress = (stress['force_lg'] / stress['area']) * 160.218  # Convert to GPa
            #stress = uniform_filter1d(stress, size=3)

            # Use unique color and linestyle for each plot
            plt.plot(
                strain, stress,
                marker='o', markersize='2',
                label=f'{key} atoms',
                color=colors[idx % len(colors)],
                linestyle=linestyles[idx % len(linestyles)],
            )
            #plt.text(strain.iloc[-1]-0.002, stress.iloc[-1]+0.02, f'{key} atoms',  color=colors[idx % len(colors)])
            idx+=1

    plt.grid(True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Strain $\epsilon$', fontsize=14, labelpad=10)
    plt.ylabel('Stress $\sigma$ (GPa)',  fontsize=14, labelpad=10)
    plt.title('Stress vs Strain',  fontsize=18)
    plt.legend()
    plt.show()
