######################### FUNCTIONS FOR FEM PDE#################################
from dolfin import *

# Differentiate
def D(v,i,z):
    if i == 2:
        return diff(v,z)
    else:
        return Dx(v,i)

# Complex number product
def c_prod(c_r,c_i,v_r,v_i,u_r,u_i):
    
    p_r  = v_r*u_r - v_i*u_i
    p_i  = v_r*u_i + v_i*u_r
    pr_r = c_r*p_r - c_i*p_i
    pr_i = c_r*p_i + c_i*p_r

    return pr_r + pr_i

# Inner product for curl of E
def v_ceps_u(vr,vi,ur,ui,ceps,z):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    S = S + c_prod(ceps(i,j,k,l),0,\
                                    D(vr[i],j,z),-D(vi[i],j,z), \
                                    D(ur[k],l,z), D(ui[k],l,z))
    return S

# Normal Inner product
def v_u(vr,vi,ur,ui):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        S = S + c_prod(1,0,vr[i],-vi[i],ur[i],ui[i])
    
    return S

# Inner product for Maxwell RHS
def v_stress_u(vr,vi,sr,si,ur,ui,P,z):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    S = S + c_prod(tensor_coeff(i,j,k,l,P),0,\
                            vr[i],-vi[i], \
                            Sij(sr,j,k,z)*ur[l] - Sij(si,j,k,z)*ui[l],\
                            Sij(sr,j,k,z)*ui[l] + Sij(si,j,k,z)*ur[l])
    return S

# Inner product for stress LHS
def v_c_ijkl_u(vr,vi,ur,ui,pr,pi,z):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    S = S + c_prod(tensor_coeff(i,j,k,l,pr),tensor_coeff(i,j,k,l,pi),\
                            D(vr[i],j,z),-D(vi[i],j,z),\
                            D(ur[k],l,z), D(ui[k],l,z))
    return S

# Inner product for stress RHS
def v_p_ijkl_EE(vr,vi,Epr,Epi,Esr,Esi,P,z):
    VEC_DIM = 3
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    S = S + c_prod(tensor_coeff(i,j,k,l,P),0,\
                            D(vr[i],j,z),-D(vi[i],j,z),\
                            Epr[k]*Esr[l] - Epi[k]*Esi[l], Epr[k]*Esr[l] - Epi[k]*Esi[l])
    return S

# contracted epsilon identity
def ceps(i,j,l,m):
    return (i==l)*(j==m) - (i==m)*(j==l)

# Sij strain tensor
def Sij(v,i,j,z):
    return 0.5*(D(v[i],j,z) + D(v[j],i,z))

# Getting Params from Tensor Cijkl/Pijkl
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

# LHS of linear elasticity equation
def v_j_sig_ji(v,u,C):
    VEC_DIM = 2
    S = 0.0
    for i in range(VEC_DIM):
        for j in range(VEC_DIM):
            for k in range(VEC_DIM):
                for l in range(VEC_DIM):
                    Sij = Dx(v[i],j)#0.5*(Dx(v[i],j) + Dx(v[j],i))
                    Skl = Dx(u[k],l)#0.5*(Dx(u[k],l) + Dx(u[l],k))
                    S = S + Sij*tensor_coeff(i,j,k,l,C)*Skl
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
