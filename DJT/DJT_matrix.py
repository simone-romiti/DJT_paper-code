# Matrix representation of the DJT

import sympy as sp
import numpy as np
from sympy.physics.quantum.spin import WignerD
import os

from S3_sphere import partition
from S3_sphere.indices import *


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
    DJT = np.matrix(np.zeros(shape=(N_alpha, N_q), dtype=complex))
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

