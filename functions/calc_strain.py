#calc_strain.py
from dolfin import *
from mesh_builder import *
from parameters import *
from matplotlib import pyplot as plt

def calc_strain(load, mesh, markers, fiber='chalco'):

    load = load

    N = 2

    # Creating the Function Space for the Strain
    ME = VectorFunctionSpace(mesh, 'Lagrange', 2)
    v  = TestFunctions(ME)
    u  = TrialFunctions(ME)

    # Defining Parameters
    e_r, rho, C, P, Qf = get_params(markers, fiber)
    if fiber=='chalco':
        print(fiber)
        E1 = 70e9 # Pa, glass
        E2 = 6.08e9 # Pa, PMMA
        v1 = 0.22 #[/]
        v2 = 0.327 #[/]
        g  = 10 #m/s^2
    else:
        print(fiber)
        E1 = 70e9 # Pa, glass
        E2 = 72.55e9 # Pa SMF clad
        v1 = 0.22 #[/]
        v2 = 0.17 # [/]
        g  = 10 #m/s^2

    r = 50e-6 # m
    L = 0.05 # m

    b = (4*load*g*r/pi/L/N * ((1 - v1**2)/E1 + (1 - v2**2)/E2))**0.5

    if b == 0:
        Pressure = 0.0
    else:
        Pressure = (2 * load *g/pi/b/L/N)/1e9

    print(b)
    print(Pressure)

    b = b/1e-6

    class Load1(UserExpression):
        def eval(self, value, x):
            if x[0]**2<b**2 - DOLFIN_EPS and x[1] < 0:
                value[0] = 0.0
                value[1] = Pressure*(1 - x[0]**2/b**2)**0.5
            else:
                value[0] = 0.0
                value[1] = 0.0
        def value_shape(self):
            return (2,)

    class Load2(UserExpression):
        def eval(self, value, x):
            if x[0]**2<b**2 - DOLFIN_EPS and x[1] > 0:
                value[0] = 0.0
                value[1] = -Pressure*(1 - x[0]**2/b**2)**0.5
            else:
                value[0] = 0.0
                value[1] = 0.0
        def value_shape(self):
            return (2,)

    # APPLY LOAD
    T1 = Load1(degree = 2)
    T2 = Load2(degree = 2)

    # LINEAR ELASTICITY
    a = v_j_sig_ji(v,u,C)*dx
    l = (T1[1]*v[1] + T2[1]*v[1])* ds # ds OVER BOUNDARY

    U = Function(ME)

    solve(a == l, U)

    return U

# LHS of linear elasticity equation
def v_j_sig_ji(v,u,C):
    VEC_DIM = 2
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    Sij = 0.5*(Dx(v[i],j) + Dx(v[j],i))
                    Skl = 0.5*(Dx(u[k],l) + Dx(u[l],k))
                    S = S + Sij*tensor_coeff(i,j,k,l,C)*Skl
    return S
