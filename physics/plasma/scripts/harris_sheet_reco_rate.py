
# export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64/:$LD_LIBRARY_PATH

from __future__ import print_function
import numpy as np
from vtk import vtkXMLPImageDataReader
from vtk.util import numpy_support
import glob
#import matplotlib
#matplotlib.use('TkAgg')
#from matplotlib import pyplot as plt

# adjust these!
path = '/home/scratch1/far/muphy2_output/whbg_10mom_256/*.pvti'
t_end = 200.
size = [50., 25.] # of the quarter
n_outputs = 100

dt = t_end/n_outputs

print("Path:", path)
pathLength = len(path)-6
files = glob.glob(path)
numFiles = len(files)

counter = 0
progress = 0.
psi = np.zeros(numFiles)
rec_rate = np.zeros(numFiles-1)

for fname in files:

	numberStr = fname[6+pathLength:]
	numberStr = numberStr[:-5]
	counter = int(numberStr)

	reader = vtkXMLPImageDataReader()
	reader.SetFileName(fname)
	reader.Update()
	img = reader.GetOutput()
	dims = [img.GetDimensions()[0]-1]

	if img.GetDimensions()[1] > 1:
	    dims.append(img.GetDimensions()[1]-1)
	if img.GetDimensions()[2] > 1:
	    dims.append(img.GetDimensions()[2]-1)

	dy = size[1]/dims[1]
	dx = dy

	# load b data
	by_data = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('By'))
	by_line = np.zeros(dims[0])

	for i in range(dims[0]):
	    by_line[i] = by_data[dims[0]*(dims[1]-1)+i]

	integral = 0.
	for i in range(0, len(by_line)):
	    integral += dx*by_line[i]

	psi[counter] = integral

	progress += 100./numFiles
	#print("Progress:", int(progress), "%", end="\r", flush=True) # python3
	print("Progress:", int(progress), "%", end="\r") # python2

for i in range(len(rec_rate)):
	rec_rate[i] = (psi[i+1]-psi[i])/dt

# output plot data
time = np.zeros(n_outputs)
for i in range(len(time)):
    time[i] = i*dt

print("time | reconnection rates")
for i in range(len(time)):
    print(time[i], rec_rate[i])

print("time | magnetic flux")
for i in range(len(time)):
    print(time[i], psi[i])

# plot
#fig = plt.figure()
#ax1 = fig.add_subplot(211)
#ax1.plot(time, rec_rate)
#ax1.set_title('Reconnection rate E_R over time')
#ax2 = fig.add_subplot(212)
#ax2.plot(time, psi[0:len(time)])
#ax2.set_title('Magnetic flux Psi over time')
#fig.subplots_adjust(hspace=.5)
#plt.show()
