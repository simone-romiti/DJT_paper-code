# Left generators in the electric basis

from su2_DJT.S3_sphere.indices import *
from su2_DJT.DJT.DJT_matrix import *

# 
def get_Lsquared(q):
    """ L^2 operator

    Args:
        q (float): half integer such that j<=q

    Returns:
        numpy matrix: \sum_a L_a*L_a = \sum_a R_a*R_a
    """
    N_q = partition.get_N_q(q)
    Lsquared = np.matrix(np.zeros(shape = (N_q, N_q)))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q)
        Lsquared[i,i] = j*(j+1)
    ####
    return Lsquared
####

# L_3
def get_L3(q):
    """L_3

    Args:
        q (float): half integer such that j<=q

    Returns:
        numpy matrix: canonical momentum L_3
    """
    N_q = partition.get_N_q(q)
    L3 = np.matrix(np.zeros(shape = (N_q, N_q)))
    for i in range(N_q):
        j, mL, mR = su2_index_to_irrep(i, q)
        L3[i,i] = mL
    ####
    return L3
####

def get_Lplus(q):
    """L_+

    Args:
        q (float): half integer such that j<=q

    Returns:
        numpy matrix: L_1 + i L_2
    """
    N_q = partition.get_N_q(q)
    Lplus = np.matrix(np.zeros(shape = (N_q, N_q)))
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
    """ L_- = (L_+)^\dagger = L_1 - i L_2 """
    return get_Lplus(q).getH()

def get_L1(q):
    """ L_1 = (L_+ + L_-)/2 """
    return (get_Lplus(q) + get_Lminus(q))/2.0
####

def get_L2(q):
    """ L_1 = (L_+ - L_-)/(2i) """
    return (get_Lplus(q) - get_Lminus(q))/(2.0*1j)
#### 

def get_La(a, q):
    """ wrapper for getting one of the L_a, a=1,2,3 """
    if a==1:
        return get_L1(q)
    elif a==2:
        return get_L2(q)
    elif a==3:
        return get_L3(q)
    else:
        raise ValueError("Invalid index of generator a={a}. Valid indices: 1,2,3".format(a=a))
####


