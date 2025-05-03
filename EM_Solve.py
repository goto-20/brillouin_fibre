from dolfin import *
from mesh_builder import *
from eigen_builder import *
from calc_strain import *
import numpy as np
import matplotlib.pyplot as plt

N_eigen = 2
path = 'E_FIELD/SMF28_4'

# Creating the mesh and obtaining the markers for core and cladding
mesh, markers = Fiber_Waveguide_Mesh(4.2,62.5)#, d_core=1+0.5,single=False)
#~ mesh, markers = Fiber_Waveguide_Mesh2()
plot(markers)
plot(mesh)
plt.show()
# Apply transverse strain
disp = calc_strain(4,mesh,markers,'silica')

e_r, _,  _,  P,  _ = get_params(markers, 'silica')

File(path + '/displacement_vector.xml') << disp

p = plot(Dx(disp[1],1) + Dx(disp[1],1))
plt.colorbar(p)
plt.savefig(path +'/eps_yy.png')
plt.show()
plt.close()

# Creating the function space for EM wave of the waveguide
combined_space, N_i, L_i, N_j, L_j = EM_Function_Generator(mesh)

# Variational Formulation
A,B, er = EM_Variational_Matrices(mesh, markers, N_i, L_i, N_j, L_j, "silica",lamb=1.55, disp = disp)

# Solving the eigenvalue for the Equation: A*u = lambda *B, where lambda= -k_z^2
E_all, n_eff = Eigenvalue_Solver(combined_space, A*dx, B*dx, er, N_eigen,dirichlet = False)


for i in range(N_eigen):
    E = E_all[i]
    plt.figure(i+1)
    p = plot(E[1])
    plt.colorbar(p)
    plt.savefig(path + '/EM_'+str(i)+'_y.png')
    File(path + '/EM_'+str(i)+'.xml') << E
    plt.close()


print(np.array([n_eff]))
np.savetxt(path + '/n.csv', np.array([n_eff]))

## Saving the mesh and markers
File(path + '/mesh.xml') << mesh
File(path + '/markers.xml') << markers

EM_Solve.py
Displaying EM_Solve.py.
