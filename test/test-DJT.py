# testing that DJT^\dagger * DJT = 1_{N_q \times N_q}

import sympy as sp
import numpy as np


import sys
sys.path.insert(0, '../')

from DJT.DJT_matrix import *

decimals = 14
print("rounding decimals", decimals)

q = 1 # truncation (integer or half integer)
print("q =", q)

N_q = partition.get_N_q(q)
N_alpha = partition.get_N_alpha(q)

N_theta = partition.get_N_theta(q)
w = partition.get_ws(N_theta)

DJT = get_DJT(q)
print(DJT.shape)
DJT_dag = get_DJT_dag(DJT = DJT)

res = DJT_dag * DJT

res = np.array(sp.Matrix(res).evalf(chop=True), dtype=complex).round(decimals=decimals)
res = np.array(res, dtype=float)
print(N_q, res.shape)

print("DJT^\dagger * DJT")
print(res)

