# Right generators in the electric basis

import partition
from DJT_matrix import *

# \sum_a R_a*R_a = \sum_a L_a*L_a
def get_Rsquared(q):
    N_q = partition.get_N_q(q)
    Rsquared = np.zeros(shape = (N_q, N_q))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q)
        Rsquared[i,i] = j*(j+1)
    ####
    return Rsquared
####

# R_3
def get_R3(q):
    N_q = partition.get_N_q(q)
    R3 = np.zeros(shape = (N_q, N_q))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q)
        R3[i,i] = mR
    ####
    return R3
####

# R_1 + i R_2
def get_Rplus(q):
    N_q = partition.get_N_q(q)
    Rplus = np.zeros(shape = (N_q, N_q))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q=q)
        if mR < j:
            i2 = su2_irrep_to_index(j, mL, mR+1, q=q)
            Rplus[i2, i] = np.sqrt(float(j*(j+1) - mR*(mR+1)))
        ####
    ####
    return Rplus
####


# R_1 - i R_2
def get_Rminus(q):
    return np.conj(get_Rplus(q)).T


def get_R1(q):
    return (get_Rplus(q) + get_Rminus(q))/(2.0)
####

def get_R2(q):
    return (get_Rplus(q) - get_Rminus(q))/(2.0*1j)
#### 

def get_Ra(a, q):
    if a==1:
        return get_R1(q)
    elif a==2:
        return get_R2(q)
    elif a==3:
        return get_R3(q)
    else:
        raise ValueError("Invalid index of generator a={a}. Valid indices: 1,2,3".format(a=a))
####




