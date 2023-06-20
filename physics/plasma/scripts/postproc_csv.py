from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from matplotlib.ticker import ScalarFormatter


print("Postprocessing of muphyII csv data")

if (len(sys.argv) > 1):
  output_directory_path = sys.argv[1]
else:
  if (sys.version_info[0] == 3):
    print("Usage:\n python3 " + sys.argv[0] + " path/to/output/directory")
  else:
    print("Usage:\n python " + sys.argv[0] + " path/to/output/directory")
  sys.exit()

path = output_directory_path + "/csv/*.csv"
data_sum = np.loadtxt(path[:-5]+"output_0_0_0.csv", delimiter=' ')
num_rows = data_sum.shape[0]
data_sum[:,1:] = 0.

files = glob.glob(path)
num_files = len(files)

i = 0
for fname in files:
  data = np.loadtxt(fname, delimiter=' ')
  num_rows = min(num_rows,data.shape[0])
  data_sum[:num_rows,1:] += data[:num_rows,1:]
  if (sys.version_info[0] > 2):
    print("Processing file", i, "out of", num_files, end="\r", flush=True) # python3
  else:
    print("Processing file", i, "out of", num_files, end="\r") # python2
  i += 1

print("\nProcessing completed")

time = data_sum[:num_rows,0]
electric_energy = data_sum[:num_rows,1]
magnetic_energy = data_sum[:num_rows,2]
em_momentum_x = data_sum[:num_rows,3]
em_momentum_y = data_sum[:num_rows,4]
em_momentum_z = data_sum[:num_rows,5]
mass_e = data_sum[:num_rows,6]
kinetic_energy_e = data_sum[:num_rows,7]
thermal_energy_e = data_sum[:num_rows,8]
momentum_x_e = data_sum[:num_rows,9]
momentum_y_e = data_sum[:num_rows,10]
momentum_z_e = data_sum[:num_rows,11]
mass_i = data_sum[:num_rows,12]
kinetic_energy_i = data_sum[:num_rows,13]
thermal_energy_i = data_sum[:num_rows,14]
momentum_x_i = data_sum[:num_rows,15]
momentum_y_i = data_sum[:num_rows,16]
momentum_z_i = data_sum[:num_rows,17]


fig, axs = plt.subplots(2, 2, figsize=(10.,5.), dpi=200, facecolor='white')

for i in range(2):
  for j in range(2):
    axs[i,j].get_yaxis().get_major_formatter().set_useOffset(False)
    axs[i,j].set_xlabel('t')

# relative mass
axs[0,0].plot(time[:], mass_i[:]/mass_i[0], color="red", linestyle="-", linewidth=3)
axs[0,0].plot(time[:], mass_e[:]/mass_e[0], color="darkblue", linestyle="-", linewidth=3)
axs[0,0].set_ylabel('Total Mass (rel)')
axs[0,0].legend(['Ions','Electrons'], loc="best", prop={'size': 10})


# relative momentum
momentum = np.sqrt((momentum_x_e[:] + momentum_x_i[:] + em_momentum_x[:])**2 + \
                   (momentum_y_e[:] + momentum_y_i[:] + em_momentum_y[:])**2 + \
                   (momentum_z_e[:] + momentum_z_i[:] + em_momentum_z[:])**2)
axs[0,1].plot(time[:], momentum[:]/momentum[0], color="black", linewidth=3)
axs[0,1].set_ylabel('Total Momentum (rel)')


# relative energy
energy = electric_energy[:] + magnetic_energy[:] + kinetic_energy_e[:] + \
          thermal_energy_e[:] + kinetic_energy_i[:] + thermal_energy_i[:]
axs[1,0].plot(time[:], energy[:]/energy[0], color="black", linewidth=3)
axs[1,0].set_ylabel('Total Energy (rel)')


# Log plot electric energy
axs[1,1].semilogy(time[:], electric_energy[:], color="black", linewidth=3)
axs[1,1].set_ylabel('Electric Energy')
# plot linear Landau damping rate
#axs[1,1].semilogy(time[:], 0.001*np.exp(-0.3066*time[:]), color="orange", linewidth=3)
 

fig.tight_layout()
plt.show()

