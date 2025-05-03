#acoustic_builder.py
from dolfin import *
import numpy as np
from parameters import *


def SBS_Variational_Matrices(mesh, markers, EM, material, n_prop, Omega, path='', lamb=1.55, s = False, disp = 0):
    # Creating the Function Space for the Acoustic Vibrations
    lagrange_order = 1
    V  = VectorElement("Lagrange", mesh.ufl_cell(),lagrange_order,dim=3)
    ME = FunctionSpace(mesh,MixedElement([V,V]))
    v  = TestFunctions(ME)
    u  = TrialFunctions(ME)

    # Defining Parameters
    e_r, rho, C, P, Qf = get_params(markers, material)

    # Permittivity in Free Space
    e_0 = Constant(8.85E-12/1E-6)

    # Momentum vector K for wavelength of 1.55 microns
    K   = Constant(4*n_prop*pi/lamb)

    A     = -rho * v_u(v,u)*dx
    B     = v_c_ijkl_u(v,u,C,0,K)*dx
    C     = v_c_ijkl_u(v,u,0,C,K)/Qf*dx
    L     = e_0 * v_T(v,P,EM,disp,K,e_r)*dx

    N_iter = len(Omega)
    Energy = np.zeros(N_iter)

    # Dirichlet condition
    zero         = Constant( (0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
    dirichlet_bc = DirichletBC(ME, zero, Zero_Boundary())

    U = Function(ME)

    # Solve for each energy E
    for i in range(N_iter):
        if i%10 == 0:
            print(path + '_' + str(100.0*i/N_iter)+'%: '+str(Omega[i])+'GHz')
        F = Constant(2*pi*Omega[i])
        M = F**2*A + B + F*C
        solve(M == L, U, dirichlet_bc)

        Ac = rho*(U[0]**2+U[1]**2+U[2]**2+U[3]**2+U[4]**2+U[5]**2)
        #~ I  = sqrt(e_r)*(EM[0]**2 + EM[1]**2 + EM[2]**2)
        #~ U_n_I = assemble(sqrt(Ac)*I*dx)
        #~ U_norm = assemble(Ac*dx)
        #~ I_norm = assemble(I**2*dx)
        #~ OVERLAP = U_n_I**2 / U_norm / I_norm

        E = assemble(F**2*Ac*dx)

        Energy[i] = E#*OVERLAP

        if s:
            File(path + '/Acoustic/'+str(Omega[i])+'_GHz.xml') << U

    return Energy

## Define Zero Boundary class
class Zero_Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Calculus on tensors :S
def tensor_coeff(i,j,k,l,param):

    if param == 0:
        return 0.0

    if (i==j) and (k==l):
        if (i==k):
            return param[0]
        else:
            return param[1]
    elif (i!=j) and (((i==k) and (j==l)) or ((i==l) and (j==k))):
        return param[2]
    else:
        return 0.0

def c_prod(c_r,c_i,v_r,v_i,u_r,u_i):

    p_r  = v_r*u_r - v_i*u_i
    p_i  = v_r*u_i + v_i*u_r
    pr_r = c_r*p_r - c_i*p_i
    pr_i = c_r*p_i + c_i*p_r

    return pr_r + pr_i

def Dr(v,K,i,j):

    if j == 2:
        return K*v[1][i]
    else:
        return Dx(v[0][i],j)

def Di(v,K,i,j):

    if j == 2:
        return -K*v[0][i]
    else:
        return Dx(v[1][i],j)

def EE_r(E,i,j):

    return E[i]*E[j]


def v_u(v,u):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        S = S + c_prod(1,0,v[0][i],-v[1][i],u[0][i],u[1][i])

    return S

def v_c_ijkl_u(v,u,pr,pi,K):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    S = S + c_prod(tensor_coeff(i,j,k,l,pr),\
                                   tensor_coeff(i,j,k,l,pi),\
                                   Dr(v,K,i,j),-Di(v,K,i,j),\
                                   Dr(u,K,k,l), Di(u,K,k,l))
    return S

def v_T(v,P,E,D,K,e_r):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    S = S + e_r**2*c_prod(tensor_coeff(k,l,i,j,P),0,\
                                   Dr(v,K,i,j),-Di(v,K,i,j),\
                                   EE_r(E,k,l),0)
                    if D != 0:
                        D_DIM = 2
                        for m in range(D_DIM):
                            for n in range(D_DIM):
                                S = S + e_r**4*eps_ij(D,i,m,P)*eps_ij(D,j,n,P)\
                                    *c_prod(tensor_coeff(k,l,m,n,P),0,\
                                            Dr(v,K,i,j),-Di(v,K,i,j),\
                                            EE_r(E,k,l),0)
                        for m in range(D_DIM):
                            S = S + e_r**3*eps_ij(D,i,m,P)\
                                *c_prod(tensor_coeff(k,l,m,j,P),0,\
                                        Dr(v,K,i,j),-Di(v,K,i,j),\
                                        EE_r(E,k,l),0)
                        for n in range(D_DIM):
                            S = S + e_r**3*eps_ij(D,j,n,P)\
                                *c_prod(tensor_coeff(k,l,i,n,P),0,\
                                        Dr(v,K,i,j),-Di(v,K,i,j),\
                                        EE_r(E,k,l),0)
    return S

# Change in Dielectric eps_ij due to strain
def eps_ij(u,i,j,P):
    VEC_DIM = 2
    S = 0.0

    for k in range(VEC_DIM):
        for l in range(VEC_DIM):
            Skl = 0.5*(Dx(u[k],l) + Dx(u[l],k))
            S = S + tensor_coeff(i,j,k,l,P)*Skl
    return -S
