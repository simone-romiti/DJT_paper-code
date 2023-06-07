# Matrix representation of the DJT

import sympy as sp
import numpy as np
from sympy.physics.quantum.spin import WignerD

try:
    from . import partition
    from .indices import *
except:
    import partition
    from indices import *


def get_DJT(q):
    N_q = partition.get_N_q(q)
    N_theta = partition.get_N_theta(q)
    N_phi = partition.get_N_phi(q)
    N_psi = partition.get_N_psi(q)
    N_alpha = partition.get_N_alpha(q)
    theta = partition.get_theta(N_theta)
    phi = partition.get_phi(N_phi)
    psi = partition.get_psi(N_psi)
    w = partition.get_ws(N_theta)
    DJT = np.zeros(shape=(N_alpha, N_q), dtype=complex)
    for i in range(N_alpha):
        i_theta, i_psi, i_phi = S3_point_to_angles_index(i, q)
        for k in range(N_q):
            j, m, mu = su2_index_to_irrep(k, q)
            s = i_theta
            D_jmmu = WignerD(sp.Rational(j), sp.Rational(m), sp.Rational(
                mu), phi[i_phi], theta[i_theta], psi[i_psi]).doit()
            DJT[i, k] = sp.Pow(j + 1/2, 1/2)*sp.sqrt(w[s]/(N_psi*N_phi))*D_jmmu
        ####
    ####
    return DJT
####

# DJT^{\dagger}
def get_DJT_dag(DJT):
    return np.conj(DJT).T
####

# def get_W(q):
#     N_theta = partition.get_N_theta(q=q)
#     N_alpha = partition.get_N_alpha(q=q)
#     w = partition.get_ws(N_theta)
#     W = np.zeros(shape=(N_alpha, N_alpha))
#     for i in range(N_alpha):
#         s = S3_point_to_angles_index(i, q)[0]  # index of theta
#         W[i, i] = np.sqrt(float(w[s]))
#     ####
#     return W
# ####

# def get_W_inv(q):
#     N_theta = partition.get_N_theta(q=q)
#     N_alpha = partition.get_N_alpha(q=q)
#     w = partition.get_ws(N_theta)
#     W_inv = np.zeros(shape=(N_alpha, N_alpha))
#     for i in range(N_alpha):
#         s = S3_point_to_angles_index(i, q)[0]  # index of theta
#         W_inv[i, i] = 1.0/np.sqrt(float(w[s]))
#     ####
#     return W_inv
# ####

# def get_V(DJT, q):
#     return np.dot(get_W_inv(q), DJT)
# ####


# def get_V_inv(DJT, q):
#     return np.dot(np.conj(DJT).T, get_W(q))
# ####

# norm squared of a numpy vector: |v|^2
# v must have shape, e.g. (100,)
def get_norm2(v):
    v_dag = np.conj(v).T
    norm2 = np.dot(v_dag, v)
    return norm2
####

# discrete version of the eigenstate of the continuum manifold (already normalized)
def get_DJT_column(DJT, j, m, mu, q):
    idx = su2_irrep_to_index(j=j, mL=m, mR=mu, q=q)
    v = DJT[:, idx]
    return v
####

# def get_WignerD_vec_values(DJT, j, m, mu, q):
#     v = get_DJT_column(j=j, m=m, mu=mu, q=q)
#     W_inv  = get_W_inv(q=q)
#     return W_inv*v
# ####

# def get_eigenstate(j, m, mu, q):
#     N_theta = partition.get_N_theta(q)
#     N_phi = partition.get_N_phi(q)
#     N_psi = partition.get_N_psi(q)
#     N_alpha = partition.get_N_alpha(q)
#     theta = partition.get_theta(N_theta)
#     phi = partition.get_phi(N_phi)
#     psi = partition.get_psi(N_psi)
#     v = np.zeros(shape=(N_alpha, 1), dtype=complex)
#     for i in range(N_alpha):
#         i_theta, i_psi, i_phi = S3_point_to_angles_index(i, q)
#         D_jmmu = WignerD(sp.Rational(j), sp.Rational(m), sp.Rational(
#             mu), phi[i_phi], theta[i_theta], psi[i_psi]).doit()
#         vi_sym = sp.Pow(j + 1/2, 1/2)*sp.sqrt(1.0/(N_psi*N_phi))*D_jmmu
#         v[i, 0] = complex(vi_sym.evalf())
#     ####
#     return v/np.sqrt(get_norm2(v))
# ####
