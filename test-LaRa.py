## testing the Lie algebra commutation relations

import sympy as sp
##from sympy.physics.matrices import msigma

import partition
from DJT_matrix import *
import operators
import electric_basis as eb

tau = operators.tau

N_c = operators.N_c
N_g = operators.N_g

q = 1
print("q =", q)

decimals = 14
print("rounding decimals", decimals)

DJT = get_DJT(q)
DJT_dag = np.conj(DJT).T

N_q = partition.get_N_q(q=q)

L = {a: operators.get_La(a=a, V=DJT, V_inv = DJT_dag, q=q) for a in [1,2,3]}
R = {a: operators.get_Ra(a=a, V=DJT, V_inv = DJT_dag, q=q) for a in [1,2,3]}
Lplus = L[1] + 1j*L[2]
Lminus = np.conj(Lplus).T
Rplus = R[1] + 1j*R[2]
Rminus = np.conj(Rplus).T

A = (np.dot(L[3], Lplus) - np.dot(Lplus, L[3])).round(decimals=decimals)
B = Lplus.round(decimals=decimals)

print(A-B)

C = (np.dot(R[3], Rplus) - np.dot(Rplus, R[3])).round(decimals=decimals)
D = Rplus.round(decimals=decimals)

print(C-D)
