
from .operators import get_U, get_La, get_Ra, get_Lsquared

from .DJT_matrix import get_V, get_V_inv, get_DJT

import numpy as np


def get_partition_list(q):
    U = get_U(q)

    partSize = U[0][0].shape[0]

    out = np.zeros((partSize, 4))

    out[:0] = U[0][0].diagonal().real()
    out[:, 1] = U[0][0].diagonal().imag()
    out[:, 2] = U[0][1].diagonal().real()
    out[:, 3] = U[0][1].diagonal().imag()

    return out


def get_canonical_momenta(q):

    DJT = get_DJT(q)

    V = get_V(DJT, q)
    V_inv = get_V_inv(DJT, q)

    L1 = get_La(1, q, V, V_inv)
    L2 = get_La(2, q, V, V_inv)
    L3 = get_La(3, q, V, V_inv)

    Lsq = get_Lsquared(q, V, V_inv)

    R1 = get_Ra(1, q, V, V_inv)
    R2 = get_Ra(2, q, V, V_inv)
    R3 = get_Ra(3, q, V, V_inv)

    return Lsq, L1, L2, L3, R1, R2, R3
