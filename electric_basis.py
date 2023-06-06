# operator representation in the truncated electric basis

import numpy as np
import sympy as sp

try:
    from .DJT_matrix import *
except:
    from DJT_matrix import *


# \sum_a L_a*L_a
def get_Lsquared(q):
    N_q = partition.get_N_q(q)
    Lsquared = np.zeros(shape = (N_q, N_q))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q)
        Lsquared[i,i] = j*(j+1)
    ####
    return Lsquared
####


def get_L3(q):
    N_q = partition.get_N_q(q)
    L3 = np.zeros(shape = (N_q, N_q))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q)
        L3[i,i] = mL
    ####
    return L3
####

# L_1 + i L2
def get_Lplus(q):
    N_q = partition.get_N_q(q)
    Lplus = np.zeros(shape = (N_q, N_q))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q=q)
        if mL < j:
            i2 = su2_irrep_to_index(j, mL+1, mR, q=q)
            Lplus[i2, i] = np.sqrt(float(j*(j+1) - mL*(mL+1)))
        ####
    ####
    return Lplus
####


# L_1 - i L2
def get_Lminus(q):
    N_q = partition.get_N_q(q)
    Lminus = np.zeros(shape = (N_q, N_q))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q=q)
        if mL > -j:        
            i2 = su2_irrep_to_index(j, mL-1, mR, q=q)
            Lminus[i2, i] = np.sqrt(float(j*(j+1) - mL*(mL-1)))
        ####
    ####
    return Lminus
####


def get_L1(q):
    N_q = partition.get_N_q(q)
    L1 = np.zeros(shape = (N_q, N_q), dtype=complex)
    Lplus = get_Lplus(q)
    Lminus = get_Lminus(q)
    L1 += (Lplus + Lminus)/(2.0)
    return L1
####

def get_L2(q):
    N_q = partition.get_N_q(q)
    L2 = np.zeros(shape = (N_q, N_q), dtype=complex)
    Lplus = get_Lplus(q)
    Lminus = get_Lminus(q)
    L2 += (Lplus - Lminus)/(2.0*1j)
    return L2
#### 

def get_La(a, q):
    if a==1:
        return get_L1(q)
    elif a==2:
        return get_L2(q)
    elif a==3:
        return get_L3(q)
    else:
        raise ValueError("Invalid index of generator a={a}. Valid indices: 1,2,3".format(a=a))
####





