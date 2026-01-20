import numpy as np
import matplotlib.pyplot as plt

# Load file and identify the last data block
with open('../output_scratch/tov.user.hst', 'r') as f:
    lines = f.readlines()

header_1 = "# Athena++ history data"
header_2 = "#  [1]=time      [2]=dt       [3]=rho-max    [4]=alpha-min"

# Find the index of the start of the last block
last_idx = 0
for i in range(len(lines) - 1):
    if header_1 in lines[i] and header_2 in lines[i+1]:
        last_idx = i + 2

# Extract columns from the last block
time, rho_max = [], []
for line in lines[last_idx:]:
    if line.startswith("#") or not line.strip(): 
        continue
    parts = line.split()
    time.append(float(parts[0]))
    rho_max.append(float(parts[2]))

# Plot the results
plt.figure(figsize=(4, 3))
plt.plot(time, np.array(rho_max) / 1e-3, color='blue', linewidth=1)
plt.xlabel(r'Time $[M_{\odot}]$', fontsize=10)
plt.ylabel(r'$\rho_{max}\ [10^{-3}M_{\odot}^{-2}]$', fontsize=10)
# plt.title('Maximum Density vs Time', fontsize=10)
# plt.legend()
plt.savefig("tov_rhomax_time.pdf", dpi=300, bbox_inches='tight')