
# export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64/:$LD_LIBRARY_PATH

from __future__ import print_function
import numpy as np
from vtk import vtkXMLPImageDataReader
from vtk.util import numpy_support
import glob
#import matplotlib
#matplotlib.use('TkAgg')
#from matplotlib import pyplot as plt


# TODO: test o-point separation for full domain


# adjust these!
path = '/home/scratch1/far/muphy2_output/isl_5mom_lambda-5_256_test/*.pvti'
lambda_ = 5.
num_outputs = 100
t_end = 2.5 # in t_A
quarter_domain = True



size = [lambda_ * np.pi, 2. * lambda_ * np.pi]
dt = size[1]*5./num_outputs

print ("Path:", path)
pathLength = len(path)-6
files = glob.glob(path)
numFiles = len(files)

counter = 0
progress = 0.
psi = np.zeros(numFiles)
island_distance = np.zeros(numFiles)
rec_rate = np.zeros(numFiles-1)
b_prime = 0.

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

	if quarter_domain == True:
	    dy = size[1]/dims[1]
	else:
	    dy = (size[1]*2.)/dims[1]
	dx = dy

	# load b data
	bx_data = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('Bx'))
	by_data = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('By'))

	if quarter_domain == True:
	    bx_line = np.zeros(dims[1])
	    by_line = np.zeros(dims[1])

	    for i in range(dims[1]):
		bx_line[i] = bx_data[dims[0]*(i+1)-1]
		by_line[i] = by_data[dims[0]*(i+1)-1]
	else:
	    bx_line = np.zeros(dims[1]/2)
	    by_line = np.zeros(dims[1]/2)

	    for i in range(dims[1]/2):
		bx_line[i] = bx_data[dims[0]*(i+1)-dims[0]/2]
		by_line[i] = by_data[dims[0]*(i+1)-dims[0]/2]

	# find o-point
	bx_line_tmp = np.abs(bx_line)
	start_idx = np.argmax(bx_line)
	end_idx = np.argmin(bx_line)
	o_point_idx = np.argmin(bx_line_tmp[start_idx:end_idx]) + start_idx
	#print(o_point_idx)
	island_distance[counter] = size[1]*(dims[1]-o_point_idx)/dims[1]*2.

	integral = 0.
	for i in range(o_point_idx, len(bx_line)):
	    integral += dy*bx_line[i]
	
	if counter == 0:
	    b_prime = np.amax(np.sqrt((bx_line[o_point_idx:])**2 + (by_line[o_point_idx:])**2))

	psi[counter] = integral

	progress += 100./numFiles
	#print ("Progress:", int(progress), "%", end="\r", flush=True) # python3
	print ("Progress:", int(progress), "%", end="\r") # python2

for i in range(len(rec_rate)):
	rec_rate[i] = (psi[i+1]-psi[i])/dt/(b_prime*b_prime)

# some output including plot data
"""
print ("b_prime:", b_prime)

print ("reconnection rates")
for i in range(len(rec_rate)):
    print (rec_rate[i])
"""

"""
print ("island distance")
for i in range(num_outputs):
    print (island_distance[i]/island_distance[0])
"""

avg_rec_rate = np.mean(rec_rate[:np.round(num_outputs*1.5/t_end)]) # average betwenn 0 and 1.5 t_A
max_rec_rate = np.amax(rec_rate)

print ("Maximum reconnection rate:", max_rec_rate)
print ("Average reconnection rate:", avg_rec_rate)

time = np.zeros(num_outputs)
for i in range(len(time)):
    time[i] = i/len(time)*t_end

"""
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(time[:np.round(num_outputs*1.5/t_end)], rec_rate[:np.round(num_outputs*1.5/t_end)])
ax1.set_title('Reconnection rate over time')
ax2 = fig.add_subplot(212)
ax2.plot(time, island_distance[:len(time)]/island_distance[0])
ax2.set_title('Distance of islands over time')
fig.subplots_adjust(hspace=.5)
plt.show()
"""
