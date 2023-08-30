# DJT could be in principle computed for each La, but it takes time when q increases

import numpy as np
import sympy as sp
from sympy.physics.quantum.spin import WignerD


import su2_DJT.operators.electric_basis as eb
import su2_DJT.S3_sphere.partition as partition
import su2_DJT.S3_sphere.indices as indices

N_c = 2  # number of colors
N_g = int(N_c*N_c - 1)  # number of generators of the Lie algebra

# Pauli matrices divided my 2
sigma_1 = np.matrix([[0, 1], [1, 0]])
sigma_2 = np.matrix([[0, -1j], [+1j, 0]])
sigma_3 = np.matrix([[1, 0], [0, -1]])

tau = {1: sigma_1/2.0, 2: sigma_2/2.0, 3: sigma_3/2.0}
tau_squared = tau[1]*tau[1] + tau[2]*tau[2] + tau[3]*tau[3]

# V* \hat{L}_a * V^{-1}, where the \hat{L}_a are the momenta in the electric basis
# V should be the DJT, and V^{-1} is DJT^{\dagger}
def get_La(a, q, V: np.matrix, V_inv: np.matrix):
    La_eb = eb.get_La(a=a, q=q)
    La = V * La_eb * V_inv
    return La
####

# \sum_a L_a*L_a. Same logic as for get_La()
def get_Lsquared(q, V: np.matrix, V_inv: np.matrix):
    Lsquared_eb = eb.get_Lsquared(q=q)
    Lsquared = V * Lsquared_eb * V_inv
    return Lsquared
    ####
####

# V* \hat{R}_a * V^{-1}, where the \hat{R}_a are the momenta in the electric basis
# V should be the DJT, and V^{-1} is DJT^{\dagger}
def get_Ra(a, q, V: np.matrix, V_inv: np.matrix):
    Ra_eb = eb.get_Ra(a=a, q=q)
    Ra = V * Ra_eb * V_inv
    return Ra
####

# \sum_a R_a*R_a. Same logic as for get_Ra()
def get_Rsquared(q, V: np.matrix, V_inv: np.matrix):
    Rsquared_eb = eb.get_Rsquared(q=q)
    Rsquared = V * Rsquared_eb * V_inv
    return Rsquared
    ####
####

# operator U in the diagonal basis.
# It has both color and Hilbert space indices,
# so we return a list of list corresponding to its matrix elements in color space
def get_U(q):
    N_theta = partition.get_N_theta(q)
    N_phi = partition.get_N_phi(q)
    N_psi = partition.get_N_psi(q)
    N_alpha = partition.get_N_alpha(q)
    theta = partition.get_theta(N_theta)
    phi = partition.get_phi(N_phi)
    psi = partition.get_psi(N_psi)
    U_11 = np.matrix(np.zeros(shape=(N_alpha, N_alpha), dtype=complex))
    U_12 = np.matrix(np.zeros(shape=(N_alpha, N_alpha), dtype=complex))
    U_21 = np.matrix(np.zeros(shape=(N_alpha, N_alpha), dtype=complex))
    U_22 = np.matrix(np.zeros(shape=(N_alpha, N_alpha), dtype=complex))
    for i in range(N_alpha):
        i_theta, i_psi, i_phi = indices.S3_point_to_angles_index(i, q)
        onehalf = sp.Rational(1/2)
        U_11[i, i] = +WignerD(onehalf, -onehalf, -onehalf,
                              phi[i_phi], theta[i_theta], psi[i_psi]).doit()
        U_12[i, i] = -WignerD(onehalf, -onehalf, +onehalf,
                              phi[i_phi], theta[i_theta], psi[i_psi]).doit()
        U_21[i, i] = -WignerD(onehalf, +onehalf, -onehalf,
                              phi[i_phi], theta[i_theta], psi[i_psi]).doit()
        U_22[i, i] = +WignerD(onehalf, +onehalf, +onehalf,
                              phi[i_phi], theta[i_theta], psi[i_psi]).doit()
    ####
    return [[U_11, U_12], [U_21, U_22]]
####

# element (a,b) in color space,
# which is a matrix of size N_alpha x N_alpha
def get_U_ab(U, a, b):
    return U[a][b]
####


# returns U^{\dagger} in the representation of get_U()
def get_Udag(U):
    # [b][a] and not [a][b] because the dagger acts also on the color space
    Udag = [[(U[b][a]).getH() for b in range(N_c)] for a in range(N_c)]
    return Udag
####


# multiplication of 2 links in color space
# if properly tensor product with the identity,
# this function can be used also for links defined at different points
def color_prod_links(U, V):
    U_prod = [[None, None], [None, None]]
    N = (U[0][0].shape)[0]
    for a in range(N_c):
        for b in range(N_c):
            U_ab = np.zeros(shape=(N, N), dtype=complex)
            for c in range(N_c):
                U_ab += (U[a][c]) * (V[c][b])
            ####
            U_prod[a][b] = U_ab
        ####
    ####
    return U_prod
####
