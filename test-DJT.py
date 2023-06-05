# testing that DJT^\dagger * DJT = 1_{N_q \times N_q}

import sympy as sp
import numpy as np


from DJT_matrix import *

decimals = 14
print("rounding decimals", decimals)

q = 1 # truncation (integer or half integer)
N_q = partition.get_N_q(q)
N_alpha = partition.get_N_alpha(q)

N_theta = partition.get_N_theta(q)
w = partition.get_ws(N_theta)

DJT = get_DJT(q)
print(DJT.shape)
DJT_H = np.conj(DJT).T

res = np.dot(DJT_H, DJT)

res = np.array(sp.Matrix(res).evalf(chop=True), dtype=complex).round(decimals=decimals)
print(N_q, res.shape)

print("DJT^\dagger * DJT")
print(res)

print("V^{-1}*V")
V = get_V(DJT = DJT, q = q)
V_inv = get_V_inv(DJT = DJT, q = q)

V_inv_V = np.dot(V_inv, V).round(decimals=decimals)
for i in range(N_q):
    for j in range(N_q):
        print(V_inv_V[i,j], end=" ")
    ####
    print("")
####