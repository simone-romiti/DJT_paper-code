# operator matrices for the spacetime lattice

import numpy as np

import su2_DJT.S3_sphere.partition as partition
from su2_DJT.lattice.grid import grid
import su2_DJT.operators.operators as operators_pt

N_c = operators_pt.N_c

# U(x, mu) operator from DJT and spacetime lattice grid of points G
# q=truncation in j
# It is implemented as a block matrix:
# [ 1  0  0   ...       0
#   0  1  0   ...       0
#   ...
#   ...    U(x, mu) ... 0
#   ...
#   0 ...            0  1 ]
# the function returns a list of 2*2=4 operators corresponding to the color indices
def get_U(G, x, mu, q):
    N_p = G.get_N() # total number of pairs (x, mu)
    N_alpha = partition.get_N_alpha(q)
    N = N_p*N_alpha
    U_11 = np.matrix(np.identity(n=N, dtype=complex))
    U_12 = np.matrix(np.identity(n=N, dtype=complex))
    U_21 = np.matrix(np.identity(n=N, dtype=complex))
    U_22 = np.matrix(np.identity(n=N, dtype=complex))
    U_pt = operators_pt.get_U(q = q)
    i = G.x_mu_to_index(x, mu)
    # print(i, q, N_p, N_alpha)
    i1, i2 = i, i + N_alpha
    U_11[i1:i2, i1:i2] = U_pt[0][0]
    U_12[i1:i2, i1:i2] = U_pt[0][1]
    U_21[i1:i2, i1:i2] = U_pt[1][0]
    U_22[i1:i2, i1:i2] = U_pt[1][1]
    ####
    U = [[U_11, U_12], [U_21, U_22]]
    return U
####

# returns U^{\dagger} in the representation of get_U()
def get_Udag(U):
    return operators_pt.get_Udag(U) # same implementation
####

def get_plaquette(G, x, mu, nu, q):
    N_p = G.get_N() # total number of pairs (x, mu)
    sizes = G.get_sizes()
    N_alpha = partition.get_N_alpha(q)
    N = N_p*N_alpha
    tr_P_munu = np.matrix(np.zeros(shape=(N, N), dtype=complex))
    U1 = get_U(G=G, x=x, mu=mu, q=q)
    x = list(x)
    x_dum = x # dummy variable
    x_dum[mu] = (x_dum[mu]+1) % sizes[mu]
    U2 = get_U(G=G, x=x_dum, mu=nu, q=q)
    x_dum = x # back to application point
    x_dum[nu] = (x_dum[nu]+1) % sizes[nu]
    U3 = get_Udag(U = get_U(G=G, x=x_dum, mu=mu, q=q))
    U4 = get_Udag(U = get_U(G=G, x=x, mu=mu, q=q))
    for a in range(N_c):
        for b in range(N_c):
            for c in range(N_c):
                for d in range(N_c):
                    tr_P_munu += (U1[a][b])*(U2[b][c])*(U3[c][d])*(U4[d][a])
                ####
            ####
        ####
    ####
    return tr_P_munu
####

# L_a(x, mu) operator from DJT and spacetime lattice grid of points G
# q=truncation in j
# It is implemented as a block matrix:
# [ 0  0  0   ...         0
#   0  0  0   ...         0
#   ...
#   ...    L_a(x, mu) ... 0
#   ...
#   0 ...              0  0 ]
def get_La(a, G, x, mu, q, V, V_inv):
    N_p = G.get_N() # total number of pairs (x, mu)
    N_alpha = partition.get_N_alpha(q)
    N = N_p*N_alpha
    La = np.matrix(np.zeros(shape=(N, N), dtype=complex))
    La_pt = operators_pt.get_La(a=q, q=q, V=V, V_inv=V_inv)
    i = G.x_mu_to_index(x, mu)
    i1, i2 = i, i + N_alpha
    La[i1:i2, i1:i2] = La_pt
    ####
    return La
####

# same as get_La() but for \sum_a L_a L_a
def get_Lsquared(G, x, mu, q, V, V_inv):
    N_p = G.get_N() # total number of pairs (x, mu)
    N_alpha = partition.get_N_alpha(q=q)
    N = N_p*N_alpha
    Lsquared = np.matrix(np.zeros(shape=(N, N), dtype=complex))
    Lsquared_pt = operators_pt.get_Lsquared(q=q, V=V, V_inv=V_inv)
    i = G.x_mu_to_index(x, mu)
    i1, i2 = i, i + N_alpha
    Lsquared[i1:i2, i1:i2] = Lsquared_pt
    ####
    return Lsquared
####
