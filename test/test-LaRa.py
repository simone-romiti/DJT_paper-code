## testing the Lie algebra commutation relations

import sympy as sp
##from sympy.physics.matrices import msigma

import sys
sys.path.insert(0, '../')

import su2_DJT.S3_sphere.partition as partition
from su2_DJT.DJT.DJT_matrix import *
import su2_DJT.operators.operators as operators
import su2_DJT.operators.electric_basis as eb

tau = operators.tau

N_c = operators.N_c
N_g = operators.N_g

q = 1
print("q =", q)

decimals = 14
print("rounding decimals", decimals)

DJT = get_DJT(q)
DJT_dag = DJT.getH()

N_q = partition.get_N_q(q=q)

L = {a: operators.get_La(a=a, V=DJT, V_inv = DJT_dag, q=q) for a in [1,2,3]}
R = {a: operators.get_Ra(a=a, V=DJT, V_inv = DJT_dag, q=q) for a in [1,2,3]}
Lplus = L[1] + 1j*L[2]
Lminus = Lplus.getH()
Rplus = R[1] + 1j*R[2]
Rminus = Rplus.getH()

A = (L[3] * Lplus - Lplus * L[3]).round(decimals=decimals)
B = Lplus.round(decimals=decimals)
d_AB = A-B

print("Check: [L_3, L_+] - L_+ = 0 ?")
print(np.array_equal(d_AB, np.zeros(shape=d_AB.shape, dtype=complex)))

C = (R[3] * Rplus - Rplus * R[3]).round(decimals=decimals)
D = Rplus.round(decimals=decimals)
d_CD = C-D

print("Check: [R_3, R_+] - R_+ = 0 ?")
print(np.array_equal(d_CD, np.zeros(shape=d_CD.shape, dtype=complex)))
