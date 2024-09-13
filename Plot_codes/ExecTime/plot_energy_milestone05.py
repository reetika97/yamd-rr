import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec

# File paths
file1 = r"run_time_berendsen.csv"
file2 = r"run_time_with_rc.csv"

# Read CSV files
runtime_b = pd.read_csv(file1, sep=';')
runtime_rc = pd.read_csv(file2, sep=';')

# Setting font to Arial
plt.rcParams["font.family"] = "Arial"

# Create a figure
fig = plt.figure(figsize=(10, 10))
gs = GridSpec(2, 2, figure=fig)

# 1. Plotting both lines together at the top (1 row, spanning both columns)
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(runtime_b['num_atoms'], runtime_b['exec_time'], color='blue', linewidth=2, linestyle='--', label='Without Rc')
ax1.plot(runtime_rc['num_atoms'], runtime_rc['exec_time'], color='red', linewidth=2, linestyle='-.', label='With Rc')
ax1.text(runtime_b['num_atoms'][20]+15, runtime_b['exec_time'][26]+15, '$Without\_R_c$', fontsize=12, color='blue')
ax1.text(runtime_rc['num_atoms'][26]-2, runtime_rc['exec_time'][26]-12, '$With\_R_c$', fontsize=12, color='red')
ax1.set_xlabel('Number of Atoms', fontsize=12)
ax1.set_ylabel('Execution time (s)', fontsize=12)
ax1.set_title('Comparative View')
ax1.grid(True)
# Add label (a)
ax1.text(0.02, 0.95, '(a)', transform=ax1.transAxes, fontsize=14, va='top', ha='left')

# 2. Plotting only the "Without Rc" line on the bottom left
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(runtime_b['num_atoms'], runtime_b['exec_time'], color='blue', linewidth=2, linestyle='--', label='Without Rc')
ax2.text(runtime_b['num_atoms'][20], runtime_b['exec_time'][26]+15, '$Without\_R_c$', fontsize=12, color='blue')
ax2.set_xlabel('Number of Atoms', fontsize=12)
ax2.set_ylabel('Execution time (s)', fontsize=12)
ax2.set_title('Without Cutoff Radius')
ax2.grid(True)
# Add label (b)
ax2.text(0.02, 0.95, '(b)', transform=ax2.transAxes, fontsize=14, va='top', ha='left')

# 3. Plotting only the "With Rc" line on the bottom right
ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(runtime_rc['num_atoms'], runtime_rc['exec_time'], color='red', linewidth=2, linestyle='-.', label='With Rc')
ax3.text(runtime_rc['num_atoms'][26]-3, runtime_rc['exec_time'][26]-15, '$With\_R_c$', fontsize=12, color='red')
ax3.set_xlabel('Number of Atoms', fontsize=12)
ax3.set_ylabel('Execution time (s)', fontsize=12)
ax3.set_title('With Cutoff Radius')
ax3.grid(True)
# Add label (c)
ax3.text(0.02, 0.95, '(c)', transform=ax3.transAxes, fontsize=14, va='top', ha='left')

# Add a main title for the entire figure
fig.suptitle('Scaling of Execution Time with Number of Atoms', fontsize=16)

# Adjust layout to prevent overlap and make room for the suptitle
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Show the plots
plt.show()