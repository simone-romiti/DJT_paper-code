# Matrix representation of the DJT

import sympy as sp
import numpy as np
from sympy.physics.quantum.spin import WignerD
import spherical_functions
import os

from su2_DJT.S3_sphere import partition
from su2_DJT.S3_sphere.indices import *


def get_WignerD(j, mL, mR, alpha, beta, gamma, use_sympy):
    if(j < 29 and not use_sympy): # spherical_functions package is stable up to this point
        # complex conjugate of the value obtained with sympy
        return np.conj(spherical_functions.Wigner_D_element(alpha, beta, gamma, j, mL, mR))
    else:
        return complex(WignerD(sp.Rational(j), sp.Rational(mL), sp.Rational(mR), alpha, beta, gamma).doit())
####

def get_DJT(q, use_sympy=False):
    N_q = partition.get_N_q(q)
    N_theta = partition.get_N_theta(q)
    N_phi = partition.get_N_phi(q)
    N_psi = partition.get_N_psi(q)
    N_alpha = partition.get_N_alpha(q)
    theta = partition.get_theta(N_theta)
    phi = partition.get_phi(N_phi)
    psi = partition.get_psi(N_psi)
    w = partition.get_ws(N_theta)
    DJT = np.matrix(np.zeros(shape=(N_alpha, N_q), dtype=complex))
    for i in range(N_alpha):
        i_theta, i_psi, i_phi = S3_point_to_angles_index(i, q)
        for k in range(N_q):
            j, m, mu = su2_index_to_irrep(k, q)
            s = i_theta
            phi_i, theta_i, psi_i = phi[i_phi], theta[i_theta], psi[i_psi]
            D_jmmu = get_WignerD(sp.Rational(j), sp.Rational(m), sp.Rational(mu), phi_i, theta_i, psi_i, use_sympy=use_sympy)
            DJT[i, k] = sp.Pow(j + 1/2, 1/2)*sp.sqrt(w[s]/(N_psi*N_phi))*D_jmmu
        ####
    ####
    return DJT
####

# DJT^{\dagger}
def get_DJT_dag(DJT):
    return DJT.getH()
####

# norm squared of a numpy vector: |v|^2
# v is a numpy.matrix --> implicit check of the dimensions when doing the product
def get_norm2(v):
    v_dag = v.getH()
    norm2 = v_dag*v
    return norm2[0,0]
####

# discrete version of the eigenstate of the continuum manifold (already normalized)
def get_DJT_column(DJT, j, m, mu, q):
    idx = su2_irrep_to_index(j=j, mL=m, mR=mu, q=q)
    v = DJT[:, idx]
    return v
####

