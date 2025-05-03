from dolfin import *
from functions.mesh_builder import *
from functions.eigen_builder import *
import numpy as np
import matplotlib.pyplot as plt
import os

path = 'E_Field'

if os.path.isdir(os.getcwd()+'/' + path) == False:
    os.mkdir(os.getcwd()+'/' + path)

N_eigen = 2
material = 'chalco'

# Creating the mesh and obtaining the markers for core and cladding
mesh, markers = Fiber_Waveguide_Mesh(1,5)
plot(markers)
plot(mesh)
plt.savefig(path + '/' + 'mesh.png')
plt.close()

# Creating the function space for EM wave of the waveguide
combined_space, N_i, L_i, N_j, L_j = EM_Function_Generator(mesh)

# Variational Formulation
A,B, er = EM_Variational_Matrices(mesh, markers, N_i, L_i, N_j, L_j, material,lamb=1.55)

# Solving the eigenvalue for the Equation: A*u = lambda *B, where lambda= -k_z^2
E_all, n_eff = Eigenvalue_Solver(combined_space, A*dx, B*dx, er, N_eigen, dirichlet = True)


for i in range(N_eigen):
    E = E_all[i]
    plt.figure(i+1)
    p = plot(E[1])
    plt.colorbar(p)
    plt.savefig(path + '/' +'EM_'+str(i)+'_y.png')
    
    plt.figure(i+2)
    p = plot(E[0])
    plt.colorbar(p)
    plt.savefig(path + '/' +'EM_'+str(i)+'_x.png')
    
    File(path + '/' +'EM_'+str(i)+'.xml') << E
    plt.close()


print(np.array([n_eff]))
np.savetxt(path + '/' +'n.csv', np.array([n_eff]))

## Saving the mesh and markers
File(path + '/' +'mesh.xml') << mesh
File(path + '/' +'markers.xml') << markers
