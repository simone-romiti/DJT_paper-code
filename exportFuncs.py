
from operators.operators import get_U, get_Lsquared, get_La  # , get_Ra

from DJT.DJT_matrix import get_V, get_V_inv, get_DJT

import numpy as np

from scipy.sparse import bmat


def get_partition_list(q):
    U = get_U(q)

    partSize = U[0][0].shape[0]

    out = np.zeros((partSize, 4))

    out[:, 0] = np.real(np.diagonal(U[0][0]))
    out[:, 1] = np.imag(np.diagonal(U[0][0]))
    out[:, 2] = np.real(np.diagonal(U[0][1]))
    out[:, 3] = np.imag(np.diagonal(U[0][1]))

    return out


def get_canonical_momenta(q):

    DJT = get_DJT(q)

    V = get_V(DJT, q)
    V_inv = get_V_inv(DJT, q)

    L1 = get_La(1, q, V, V_inv)
    L2 = get_La(2, q, V, V_inv)
    L3 = get_La(3, q, V, V_inv)

    Lsq = get_Lsquared(q, V, V_inv)

    R1 = get_La(1, q, V, V_inv)
    R2 = get_La(2, q, V, V_inv)
    R3 = get_La(3, q, V, V_inv)


    # out = [bmat([[Lsq]], format="csc"), bmat([[L1]], format="csc"), bmat([[L2]], format="csc"), bmat([[L3]], format="csc"), bmat([[R1]], format="csc"), bmat([[R2]], format="csc"), bmat([[R3]], format="csc")]

    out = [Lsq, L1, L2, L3, R1, R2, R3]

    return out
