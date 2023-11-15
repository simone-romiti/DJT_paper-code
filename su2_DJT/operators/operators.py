"""
Operators wrappers

Note: DJT could be in principle computed for each La, but it takes time when q increases
"""


import numpy as np
import sympy as sp

import su2_DJT.operators.electric_basis as eb
import su2_DJT.S3_sphere.partition as partition
import su2_DJT.S3_sphere.indices as indices
from su2_DJT.DJT.DJT_matrix import get_WignerD

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
    """
    V* \hat{L}_a * V^{-1}, where the \hat{L}_a are the momenta in the electric basis
    V is the DJT, and V^{-1} is DJT^{\dagger}

    Args:
        a (int): index of canonical momentum
        q (float): j<=q
        V (np.matrix): DJT matrix
        V_inv (np.matrix): DJT^\dagger

    Returns:
        np.matrix: representation of the L_a in the magnetic basis
    """
    La_eb = eb.get_La(a=a, q=q)
    La = V * La_eb * V_inv
    return La
####

def get_Lsquared(q, V: np.matrix, V_inv: np.matrix):
    """ \sum_a L_a*L_a. Same logic as for get_La() """
    Lsquared_eb = eb.get_Lsquared(q=q)
    Lsquared = V * Lsquared_eb * V_inv
    return Lsquared
    ####
####

def get_Ra(a, q, V: np.matrix, V_inv: np.matrix):
    """
    V* \hat{R}_a * V^{-1}, where the \hat{R}_a are the momenta in the electric basis
    V is the DJT, and V^{-1} is DJT^{\dagger}

    Args:
        a (int): index of canonical momentum
        q (float): j<=q
        V (np.matrix): DJT matrix
        V_inv (np.matrix): DJT^\dagger

    Returns:
        np.matrix: representation of the R_a in the magnetic basis
    """
    Ra_eb = eb.get_Ra(a=a, q=q)
    Ra = V * Ra_eb * V_inv
    return Ra
####

def get_Rsquared(q, V: np.matrix, V_inv: np.matrix):
    """ \sum_a R_a*R_a. Same logic as for get_Ra() """
    Rsquared_eb = eb.get_Rsquared(q=q)
    Rsquared = V * Rsquared_eb * V_inv
    return Rsquared
    ####
####

def get_U(q):
    """
    operator U in the diagonal basis (i.e. magnetic basis).
    It has both color and Hilbert space indices,

    Args:
        q (float): maximum value of j defining N_\alpha involved in the DJT construction

    Returns:
        list of list : each element corresponds a pair of indices in color space
    """
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
        U_11[i, i] = +get_WignerD(onehalf, -onehalf, -onehalf, phi[i_phi], theta[i_theta], psi[i_psi], use_sympy=False)
        U_12[i, i] = -get_WignerD(onehalf, -onehalf, +onehalf, phi[i_phi], theta[i_theta], psi[i_psi], use_sympy=False)
        U_21[i, i] = -get_WignerD(onehalf, +onehalf, -onehalf, phi[i_phi], theta[i_theta], psi[i_psi], use_sympy=False)
        U_22[i, i] = +get_WignerD(onehalf, +onehalf, +onehalf, phi[i_phi], theta[i_theta], psi[i_psi], use_sympy=False)
    ####
    return [[U_11, U_12], [U_21, U_22]]
####

def get_U_ab(U, a, b):
    """ Element (a,b) of the operator U, in color space.
        Note: for each pair (a,b) we get a matrix of size N_alpha x N_alpha
    """
    return U[a][b]
####


def get_Udag(U):
    """ returns U^{\dagger} in the representation of get_U() """
    # [b][a] and not [a][b] because the dagger acts also on the color space
    Udag = [[(U[b][a]).getH() for b in range(N_c)] for a in range(N_c)]
    return Udag
####


def color_prod_links(U, V):
    """
    multiplication of 2 links in color space
    if properly tensor product with the identity,
    this function can be used also for links defined at different points
    """
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
