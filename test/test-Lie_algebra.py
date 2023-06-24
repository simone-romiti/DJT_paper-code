## testing the Lie algebra commutation relations

import sympy as sp
##from sympy.physics.matrices import msigma

import sys
sys.path.insert(0, '../')

import S3_sphere.partition as partition
from DJT.DJT_matrix import *
import operators.operators as operators
import operators.electric_basis as eb

tau = operators.tau

N_c = operators.N_c
N_g = operators.N_g

q = 1
print("q =", q)

decimals = 14
print("rounding decimals", decimals)

DJT = get_DJT(q)
DJT_dag = get_DJT_dag(DJT = DJT)

N_q = partition.get_N_q(q=q)

eb_L1 = eb.get_La(a=1, q=q) 
eb_L2 = eb.get_La(a=2, q=q) 
eb_L3 = eb.get_La(a=3, q=q)
eb_L = {1: eb_L1, 2: eb_L2, 3: eb_L3}
eb_Lplus, eb_Lminus = eb.get_Lplus(q=q), eb.get_Lminus(q=q)

L1 = operators.get_La(a=1, V=DJT, V_inv = DJT_dag, q=q) 
L2 = operators.get_La(a=2, V=DJT, V_inv = DJT_dag, q=q) 
L3 = operators.get_La(a=3, V=DJT, V_inv = DJT_dag, q=q)
L = {1: L1, 2: L2, 3: L3}


eb_R1 = eb.get_Ra(a=1, q=q) 
eb_R2 = eb.get_Ra(a=2, q=q) 
eb_R3 = eb.get_Ra(a=3, q=q)
eb_R = {1: eb_R1, 2: eb_R2, 3: eb_R3}
eb_Rplus, eb_Rminus = eb.get_Lplus(q=q), eb.get_Lminus(q=q)

R1 = operators.get_Ra(a=1, V=DJT, V_inv = DJT_dag, q=q) 
R2 = operators.get_Ra(a=2, V=DJT, V_inv = DJT_dag, q=q) 
R3 = operators.get_Ra(a=3, V=DJT, V_inv = DJT_dag, q=q)
R = {1: R1, 2: R2, 3: R3}


print("Is \hat{L}_i = V^{-1} * L_i * V ?")
for i in [x+1 for x in range(N_g)]:
    delta = eb_L[i] - DJT_dag * L[i] * DJT
    delta = delta.round(decimals=decimals)
    bg = np.array_equal(delta, np.zeros(shape=(N_q, N_q), dtype=complex))
    print("i :", i, bg)
####


for i in [x+1 for x in range(N_g)]:
    for j in [x+1 for x in range(N_g)]:
        for k in [x+1 for x in range(N_g)]:
            if i==j or i==k or j==k:
                continue
            ####
            L_list = [eb_L, L]
            basis_list = ["electric", "magnetic"]
            for ib in range(2):
                A = L_list[ib][i] * L_list[ib][j] - L_list[ib][j] * L_list[ib][i]
                B = 1j*float(sp.LeviCivita(i,j,k))*L_list[ib][k]
                N = B.shape[0]
                dAB = np.array(sp.Matrix(A-B).evalf(decimals, chop=True), dtype=complex).round(decimals=decimals)
                dAB = np.array(dAB) ## shoud be the 0 matrix
                print(basis_list[ib], "basis.", "i, j, k :", i,j,k, end = "  ")
                b1 = np.array_equal(dAB, np.zeros(shape=(N, N), dtype=complex))
                if not b1:
                    print(dAB)
                    quit()
                print("[L_i,L_j] = i \espsilon_{ijk} L_k ?", b1)
            ####
            x = tau[i] * tau[j] - tau[j] * tau[i]
            y = 1j*float(sp.LeviCivita(i,j,k))*tau[k]
            dxy = np.array(sp.Matrix(x-y).evalf(chop=True), dtype=complex).round(decimals=decimals)
            dxy = np.array(dxy) ## shoud be the 0 matrix
            b2 = np.array_equal(dxy, np.zeros(shape=(N_c, N_c), dtype=complex))
            print("i, j, k :", i,j,k, end = " ")
            print("[\\tau_i,\\tau_j] = i \espsilon_{ijk} \\tau_k ?        ", b2)
        ####
    ####
####