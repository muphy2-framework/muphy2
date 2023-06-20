from __future__ import print_function
import numpy as np
from vtk import vtkXMLPImageDataReader
from vtk.util import numpy_support
import glob
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import rc


# adjust these
loadpath = 'ot_10mom_2048x_gradT/output50.pvti'
c0 = 18.174 # 20.


print("Path:", loadpath)
img_file_name = 'postproc_divEB.png'
eps0 = 1./(c0*c0)

rc('text', usetex=True)
rc('font', family='serif')

fig,axs = plt.subplots(1, 2, sharey=True, figsize=(10.,3.5))
for i in range(2):
  axs[i].set_xlabel(r'$x / d_{i,0}$')
axs[0].set_ylabel(r'$y / d_{i,0}$')
axs[0].set_title(r'$\nabla \cdot \mathbf{B}$')
axs[1].set_title(r'$\nabla \cdot \mathbf{E} - \rho / \epsilon_0$')

reader = vtkXMLPImageDataReader()
reader.SetFileName(loadpath)
reader.Update()
img = reader.GetOutput()

dims = [img.GetDimensions()[0]-1,img.GetDimensions()[1]-1,img.GetDimensions()[2]-1]
origin = img.GetOrigin()
spacing = img.GetSpacing()

dx = spacing[0]
dy = spacing[1]
x = np.array([origin[0] + i*dx for i in range(1,dims[0]-1)])
y = np.array([origin[1] + i*dy for i in range(1,dims[1]-1)])

# load fields
Ex = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('Ex'))
Ex = np.reshape(Ex,(dims[1],dims[0]))
Ey = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('Ey'))
Ey = np.reshape(Ey,(dims[1],dims[0]))
Ez = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('Ez'))
Ez = np.reshape(Ez,(dims[1],dims[0]))
Bx = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('Bx'))
Bx = np.reshape(Bx,(dims[1],dims[0]))
By = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('By'))
By = np.reshape(By,(dims[1],dims[0]))
Bz = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('Bz'))
Bz = np.reshape(Bz,(dims[1],dims[0]))
n_e = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('n_e'))
n_e = np.reshape(n_e,(dims[1],dims[0]))
n_i = numpy_support.vtk_to_numpy(img.GetCellData().GetArray('n_i'))
n_i = np.reshape(n_i,(dims[1],dims[0]))

# div B
divB = np.zeros([dims[1],dims[0]])
divB = (Bx-np.roll(Bx,1,axis=1))/dx + (By-np.roll(By,1,axis=0))/dy
divB = divB[1:-1,1:-1]

z_min = -np.maximum(np.abs(divB.min()), np.abs(divB.max()))
z_max = -z_min
plt.axis([x.min(), x.max(), y.min(), y.max()])
im = axs[0].pcolor(x, y, divB, cmap='coolwarm', vmin=z_min, vmax=z_max, rasterized=True, shading='auto')
fig.colorbar(im, ax=axs[0], pad = 0.06, aspect=30)#, shrink = .7)

# div E - rho/espilon_0
divE_minus_rho = np.zeros([dims[1],dims[0]])
divE_minus_rho = (np.roll(Ex,-1,axis=1)-Ex)/dx + (np.roll(Ey,-1,axis=0)-Ey)/dy - (n_i-n_e)/eps0
divE_minus_rho = divE_minus_rho[1:-1,1:-1]

z_min = -np.maximum(np.abs(divE_minus_rho.min()), np.abs(divE_minus_rho.max()))
z_max = -z_min
plt.axis([x.min(), x.max(), y.min(), y.max()])
im = axs[1].pcolor(x, y, divE_minus_rho, cmap='coolwarm', vmin=z_min, vmax=z_max, rasterized=True, shading='auto')
fig.colorbar(im, ax=axs[1], pad = 0.06, aspect=30)#, shrink = .7)


print("max(abs(div(E)-rho/eps_0) =", np.max(np.abs(divE_minus_rho)))
print("max(abs(div(B))) =", np.max(np.abs(divB)))


plt.savefig(img_file_name, bbox_inches='tight', dpi=300)
#plt.show()

