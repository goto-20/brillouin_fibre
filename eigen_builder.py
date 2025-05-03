#eigen_builder.py

from dolfin import *
from petsc4py import PETSc
from parameters import *
from tensor_funcs import *

def EM_Function_Generator(mesh):

    # Define function space using FEnICS
    nedelec_order = 1
    lagrange_order = 1
    V_N = FiniteElement("Nedelec 1st kind H(curl)", mesh.ufl_cell(), nedelec_order)
    V_L = FiniteElement("Lagrange", mesh.ufl_cell(), lagrange_order)
    combined_space = FunctionSpace(mesh,V_N*V_L)

    # Define basis and bilinear form
    (N_i,L_i) = TestFunctions(combined_space)
    (N_j,L_j) = TrialFunctions(combined_space)

    return combined_space, N_i, L_i, N_j, L_j

def curl_t(w):
    return Dx(w[1],0) - Dx(w[0],1)

def eps_ij_EiEj(N_v, N_u, disp, P, e_r):
    D_DIM = 2
    S = 0
    for i in range(D_DIM):
        for j in range(D_DIM):
            S = S + e_r**2*eps_ij(disp,i,j,P)*N_v[i]*N_u[j]
    return S

def delta(i,j):
    if i == j:
        return 1
    else:
        return 0

def EM_Variational_Matrices(mesh, marker, N_i, L_i, N_j, L_j, material, lamb=1.55, disp = 0):

    '''
    Setting up the variational formulation for Maxwell's equation on a 2D
    cross-section for an arbitrary waveguide
    '''

    e_r, _,  _,  P,  _ = get_params(marker, material)

    # Define Wavevector for light of lambda = 1.55 microns
    k_o_squared = Constant((2*pi/lamb)**2)

    # Relative permiability of Chalc_PMMA (=1)
    one_over_u_r = Constant(1.0)

    s_tt = one_over_u_r*dot(curl_t(N_i), curl_t(N_j))
    if disp != 0 :
        t_tt = e_r*dot(N_i, N_j) + eps_ij_EiEj(N_i, N_j, disp, P, e_r)
    else:
        t_tt = e_r*dot(N_i, N_j)
    s_zz = one_over_u_r*dot(grad(L_i), grad(L_j))
    t_zz = e_r*L_i*L_j

    b_tt = one_over_u_r*dot(N_i, N_j)
    b_tz = one_over_u_r*dot(N_i, grad(L_j))
    b_zt = one_over_u_r*dot(grad(L_i), N_j)
    b_zz = ( s_zz - k_o_squared * t_zz)

    a = ( s_tt - k_o_squared * t_tt )
    b = ( b_tt + b_tz + b_zt + b_zz )

    return a,b, e_r

class Zero_Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

def Eigenvalue_Solver(combined_space, a, b, e_r, n_eigen, dirichlet = True):

    '''
    Solving the eigenvalue problem A*u = lambda B using PETSc
    Save the eigenvector into a .pvd file
    '''

    # Test for PETSc and SLEPc
    if not has_linear_algebra_backend("PETSc"):
        print("DOLFIN has not been configured with PETSc. Exiting.")
        exit()

    if not has_slepc():
        print("DOLFIN has not been configured with SLEPc. Exiting.")
        exit()

    # Assembling the equations as Matrices
    A = PETScMatrix()
    B = PETScMatrix()
    assemble(a, tensor=A)
    assemble(b, tensor=B)

    # Setting up dirichlet Boundary
    if dirichlet:
        zero = Constant( (0.0, 0.0, 0.0) )
        dirichlet_bc = DirichletBC(combined_space, zero, Zero_Boundary())
        dirichlet_bc.apply(A)
        dirichlet_bc.apply(B)

    # Create eigensolver with generalized Hermitian condition
    # Documentation of Eigensolver : pydoc dolfin.SLEPcEigenSolver
    eigensolver = SLEPcEigenSolver(A,B)
    eigensolver.parameters["spectrum"] = "smallest real"
    #eigensolver.parameters["solver"] = "lanczos"

    # Compute all eigenvalues of A = \lambda B_RHS x
    print("Computing eigenvalues. This can take a minute.")
    eigensolver.solve(4)

    fields = []
    n_eff = []
    for i in range(n_eigen):
        try:
            (lr, lc) = eigensolver.get_eigenvalue(i)
        except RuntimeError:
            print(' Aborting Mission ')
            exit()
        print("Cutoff frequency:", lr, "n_eff:", sqrt(-lr)*1.55/2/pi)
        ## Extract eigenvector for cutoff frequency
        (_, _, rx, _) = eigensolver.get_eigenpair(i)

        #~ ur = Function(combined_space,rx)

        ## Convert from e_t to E_t by diving each term by 1/sqrt(-lr)
        ## Also neglecting axial component
        #~ U = ur.vector()
        dofs_is = [combined_space.sub(i).dofmap().dofs() for i in range(2)]
        rx[dofs_is[0][:]] *= 1/sqrt(-lr)
        rx[dofs_is[1][:]] *= 1

        ## Initialize function and assign eigenvector
        ur = Function(combined_space, rx)

        ## normalize field to unit Power
        factor = assemble(sqrt(e_r)*(ur[0]**2 + ur[1]**2 + ur[2]**2)*dx)
        ur     = Function(combined_space, rx/sqrt(factor))
        fields.append(ur)
        n_eff.append(sqrt(-lr)*1.55/2/pi)

    return fields, n_eff
