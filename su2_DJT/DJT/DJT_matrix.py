# Matrix representation of the DJT

import sympy as sp
import numpy as np
from sympy.physics.quantum.spin import WignerD
import spherical_functions
import os

from su2_DJT.S3_sphere import partition
from su2_DJT.S3_sphere.indices import *


def get_WignerD(j, mL, mR, alpha, beta, gamma, use_sympy):
    """
    Returns the Wigner D function:

    D^{j}_{mL,mR} = e^{+i*alpha} * d^{j}_{mL,mR}(beta) * e^{+i*gamma}
    
    Note: the factor (-1)^(mL-mR) is needed to compensate the different normalization 

    """
    if(j < 29 and not use_sympy): # spherical_functions package is stable up to this point
        # complex conjugate to get the same as sympy
        return (-1)**(mL-mR) * np.conj(spherical_functions.Wigner_D_element(alpha, beta, gamma, j, -mL, -mR))
    else:
        return (-1)**(mL-mR) * complex(WignerD(sp.Rational(j), -sp.Rational(mL), -sp.Rational(mR), alpha, beta, gamma).doit())
####

def get_DJT(q, use_sympy=False):
    """Returns the matrix representation of the DJT

    Args:
        q (float): half-integer such that j<=q
        use_sympy (bool, optional): Use sympy to compute the Wigner D functions. Defaults to False.

    Returns:
        numpy.matrix: matrix representation of the transform
    """
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
    """Returns the dagger of the DJT

    Args:
        DJT (numpy matrix): DJT matrix representation

    Returns:
        numpy matrix: DJT^\dagger
    """
    return DJT.getH()
####

# norm squared of a numpy vector: |v|^2
# v is a numpy.matrix --> implicit check of the dimensions when doing the product
def get_norm2(v):
    """Norm squared of a vector

    Args:
        v (numpy (column) matrix): vector representing the wavefunction at the sampling points

    Returns:
        float: ||v||^2
    """
    v_dag = v.getH()
    norm2 = v_dag*v
    return norm2[0,0]
####

def get_norm(v):
    """Norm of the vector

    Args:
        v (numpy matrix): <see get_norm2()>

    Returns:
        float: ||v||
    """
    return np.sqrt(get_norm2(v))
####

def get_DJT_column(DJT, j, m, mu, q):
    """Column of the DJT matrix corresponding to the su(2) irrep (j, m, \mu)
    Physically this is a discrete version of the eigenstate of the continuum manifold (already normalized)

    Args:
        DJT (_type_): _description_
        j (float): half integer j<=q
        m (float): half integer between -j and j
        mu (float): half integer between -j and j
        q (float): maximum value of j defining the DJT

    Returns:
        numpy matrix: DJT^i_{j,m,\mu}, i=1,\ldots,N_\alpha  
    """
    idx = su2_irrep_to_index(j=j, mL=m, mR=mu, q=q)
    v = DJT[:, idx]
    return v
####

def get_1mP_garbage_space(DJT, q):
    """ projector (1-P), where P projects to the garbage space with j>q """
    N_alpha = partition.get_N_alpha(q)
    N_q = partition.get_N_q(q)
    P = np.matrix(np.zeros(shape=(N_alpha, N_alpha), dtype=complex))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q=q)
        vj = get_DJT_column(DJT=DJT, j=j, m=mL, mu=mR, q=q)
        P += vj*(vj.getH())
    ####
    return P
####

def get_P_garbage_space(DJT, q):
    """ projector P to the garbage space with j>q """
    N_alpha = partition.get_N_alpha(q)
    Id = np.matrix(np.eye(N_alpha))
    return Id - get_1mP_garbage_space(DJT, q)
####

