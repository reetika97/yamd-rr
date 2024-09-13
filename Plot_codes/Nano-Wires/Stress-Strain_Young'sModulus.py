from matplotlib import pyplot as plt
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

if __name__ == '__main__':
    file = r'strain_stress_files/3050_0.3_0.0_svf.csv'
    plt.rcParams["font.family"] = "Arial"
    plt.figure()

    stress = pd.read_csv(file, sep=';')
    strain = stress['strain']
    #strain = uniform_filter1d(strain, size=1)
    stress = (stress['force_lg'] / stress['area']) * 160.218  # Convert to GPa
    #stress = uniform_filter1d(stress, size=1)
    #plt.plot([strain[1], strain[1], strain[0]], [stress[1], stress[0], stress[0]], color='purple', linestyle='--')
    #plt.fill([strain[1], strain[1], strain[0]], [stress[1], stress[0], stress[0]], facecolor="white", edgecolor="purple", linestyle='--', hatch="//")
    plt.plot(strain, stress, '-o', color='darkblue', linestyle='-')
    Youngs_mod=(stress[2] - stress[0]) / (strain[2] - strain[0])
    plt.text(strain[2], stress[0], f'Young\'s Modulus = {Youngs_mod:.2f}GPa', color='Purple')
    plt.grid(True)
    plt.xlabel(' strain ')
    plt.ylabel('Stress (GPa)')
    plt.title('Stress vs Strain')
    plt.show()
