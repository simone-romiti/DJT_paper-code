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

q = 3/2
print("q =", q)

decimals = 14
print("rounding decimals", decimals)

DJT = get_DJT(q)
DJT_dag = np.conj(DJT).T
V = get_V(DJT=DJT, q=q)
V_inv = get_V_inv(DJT=DJT, q=q)

N_q = partition.get_N_q(q=q)

eb_L1 = eb.get_La(a=1, q=q) 
eb_L2 = eb.get_La(a=2, q=q) 
eb_L3 = eb.get_La(a=3, q=q)
eb_L = {1: eb_L1, 2: eb_L2, 3: eb_L3}
eb_Lplus, eb_Lminus = eb.get_Lplus(q=q), eb.get_Lminus(q=q)

L1 = operators.get_La(a=1, V=V, V_inv = V_inv, q=q) 
L2 = operators.get_La(a=2, V=V, V_inv = V_inv, q=q) 
L3 = operators.get_La(a=3, V=V, V_inv = V_inv, q=q)
L = {1: L1, 2: L2, 3: L3}

print("Is \hat{L}_i = V^{-1} * L_i * V ?")
for i in [x+1 for x in range(N_g)]:
    delta = eb_L[i] - np.dot(V_inv, np.dot(L[i], V))
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
                A = np.dot(L_list[ib][i], L_list[ib][j]) - np.dot(L_list[ib][j], L_list[ib][i])
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
            x = np.dot(tau[i], tau[j]) - np.dot(tau[j], tau[i])
            y = 1j*float(sp.LeviCivita(i,j,k))*tau[k]
            dxy = np.array(sp.Matrix(x-y).evalf(chop=True), dtype=complex).round(decimals=decimals)
            dxy = np.array(dxy) ## shoud be the 0 matrix
            b2 = np.array_equal(dxy, np.zeros(shape=(N_c, N_c), dtype=complex))
            print("i, j, k :", i,j,k, end = " ")
            print("[\\tau_i,\\tau_j] = i \espsilon_{ijk} \\tau_k ?        ", b2)
        ####
    ####
####