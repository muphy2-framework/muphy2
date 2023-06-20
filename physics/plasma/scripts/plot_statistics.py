from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import glob

path = "/home/scratch1/far/muphy2_output/test/csv/*.csv"

data_sum = np.loadtxt(path[:-5]+"output_0_0_0.csv", delimiter=' ')
num_rows = data_sum.shape[0]
num_columns = data_sum.shape[1]
data_sum[:,1:] = 0.

files = glob.glob(path)
num_files = len(files)

i = 0
for fname in files:
  data = np.loadtxt(fname, delimiter=' ')
  data_sum[:,1:] += data[:,1:]
  print("Processing file", i, "out of", num_files, end="\r") # python2
  #print("Processing file", i, "out of", num_files, end="\r", flush=True) # python3
  i += 1

print("\nProcessing completed")

plt.semilogy(data_sum[:,0], data_sum[:,1])
#plt.plot(data_sum[:,0], data_sum[:,9])

plt.xlabel('t')
plt.ylabel('E')
plt.show()
