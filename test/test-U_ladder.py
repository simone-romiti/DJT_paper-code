# test the ladder action of U on the vacuum

import sympy as sp

import operators.operators as operators
from DJT.DJT_matrix import *

import operators

N_c = operators.N_c
N_g = operators.N_g

q = 1
print("q =", q)
DJT = get_DJT(q)
DJT_dag = get_DJT_dag(DJT=DJT)

U = operators.get_U(q = q)
U_dag = operators.get_Udag(U)

La = [operators.get_La(a=a, q=q, V=DJT, V_inv=DJT_dag) for a in [1, 2, 3]] 
Lsquared = operators.get_Lsquared(q=q, V=DJT, V_inv=DJT_dag)

Ra = [operators.get_Ra(a=a, q=q, V=DJT, V_inv=DJT_dag) for a in [1, 2, 3]] 
N_alpha = partition.get_N_alpha(q)

decimals = 15
print("rounding decimals", decimals)

vac = get_DJT_column(DJT, 0, 0, 0, q=q)

for jd in range(1, int(2*q)+1):
    j = jd/2
    print("j =", sp.Rational(j))
    for a in range(N_c):
        for b in range(N_c):
            v = U[a][b] * vac
            for k in range(jd-1):
                v = U[a][b] * v
            ####
            w = Lsquared * v
            dw = w - j*(j+1)*v
            msg = "| L^2 (U_{a}{b})^{j}|0> - j(j+1) (U_{a}{b})^{j}|0> |^2 =".format(j=jd, a=a, b=b)
            res = get_norm2(dw).round(decimals = decimals)
            print(msg, res)
        ####
    ####
####

