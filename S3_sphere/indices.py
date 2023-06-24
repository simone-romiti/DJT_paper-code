
import sympy as sp
import numpy as np

from . import partition


def check_su2_irrep(j, mL, mR):
    if abs(mL) > j or abs(mR) > j:
        err_mess = "Error: (j, mL, mR)=({j}, {mL}, {mR}) is not a valid su(2) irrep".format(j=j, mL=mL, mR=mR)
        raise ValueError(err_mess)
    ####
####

# (j, mL, mR) from the checkerboard index of the irrep with truncation q in j
def su2_index_to_irrep(k, q):
    i = 0
    for j1 in [q - j_i/2 for j_i in range(0, int(2*q) + 1)]:
        deg_j1 = int(2*j1) + 1 # degeneracy of the left and right quantum numbers
        for m in [j1 - i_m for i_m in range(deg_j1)]:
            for mu in [j1 - i_mu for i_mu in range(deg_j1)]:
                if i==k:
                    return [sp.Rational(j1), sp.Rational(m) , sp.Rational(mu)] 
                else:
                    i += 1
                ####
            ####
        ####
    ####
    err_mess = "{k} does not correspond to any (j, m, mu) su(2) irrep with truncation q = {q}".format(k=k, q=q)
    raise LookupError(err_mess)
####

# checkerboard index of the irrep, going from 0,...,N_q-1
# Each block of given "j" has size equal to get_N_q(j)
# Inside each block, mL and mR are the row and column index of a table of size (2*j + 1)
# Out convention is that we order the irreps from the highest to the lowest
def su2_irrep_to_index(j, mL, mR, q):
    check_su2_irrep(j, mL, mR)
    N_q = partition.get_N_q(q)
    k0 = partition.get_N_q(j)
    i0 = N_q - k0
    deg_j = int(2*j) + 1
    i1 = deg_j*(j - mL) + (j - mR)
    return int(i0+i1)
####

# i = i_theta*N_phi*N_psi + i_phi*N_psi + i_psi
def S3_point_to_angles_index(i, q):
    N_phi =   partition.get_N_phi(q)
    N_psi =   partition.get_N_psi(q)
    i_theta = int(i/(N_phi*N_psi))
    i1 = i - i_theta*N_phi*N_psi
    i_phi = int(i1/N_psi)
    i2 = i1 - i_phi*N_phi
    i_psi = i2
    return [i_theta, i_phi, i_psi]
####


# i = i_theta*N_phi*N_psi + i_psi*N_phi + i_phi
def angles_to_S3_point_index(i_theta, i_phi, i_psi, q):
    N_phi =   partition.get_N_phi(q)
    N_psi =   partition.get_N_psi(q)
    i = i_theta*N_phi*N_psi + i_phi*N_psi + i_psi
    return i
####
