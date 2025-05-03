from functions.acoustic_builder import *
from functions.eigen_builder import EM_Function_Generator
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os

folder_A = 'A_Field'
folder_E = 'E_Field'

if os.path.isdir(os.getcwd()+'/' + folder_A) == False:
    os.mkdir(os.getcwd()+'/' + folder_A)

# Import mesh and marker
mesh = Mesh(folder_E + '/mesh.xml')
marker = MeshFunction('size_t', mesh, folder_E + '/markers.xml')
n_eff = np.loadtxt(folder_E + '/n.csv')
n_eff = float(n_eff[0])
# Define Function space from mesh
combined_space = EM_Function_Generator(mesh)
E1             = Function(combined_space[0], folder_E + '/EM_0.xml')

# Frequency Axis to compute the SBS profile
N_iter = 51
Omega  = np.linspace(7.5,7.8,N_iter) # 51 points from 7.5 to 7.8 GHz

Energy = SBS_Variational_Matrices(mesh, marker, E1, "chalco", n_eff , Omega, folder_A, 1.55)

L_om = np.argmax(Energy)
print(Omega[L_om])
plt.plot(Omega, Energy, 'b')#, Omega, L, 'g')
data = [[Omega[i],Energy[i]] for i in range(len(Omega))]
np.savetxt(folder_A + '/Energy.csv', data, delimiter=',')
plt.savefig(folder_A + '/Energy.png')

plt.show()
plt.close()
