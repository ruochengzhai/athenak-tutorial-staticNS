import numpy as np
import matplotlib.pyplot as plt

# Load file and identify the last data block
data = np.genfromtxt('path/to/tov.user.hst', comments="#", usecols=(0,2)).T
time, rho_max = data[0,:], data[1,:]

# Plot the results
plt.figure(figsize=(4, 3))
plt.plot(time, rho_max / 1e-3, color='blue', linewidth=1)
plt.xlabel(r'Time $[M_{\odot}]$', fontsize=10)
plt.ylabel(r'$\rho_{max}\ [10^{-3}M_{\odot}^{-2}]$', fontsize=10)
plt.grid(True)
plt.savefig("tov_rhomax_time.pdf", dpi=300, bbox_inches='tight')
