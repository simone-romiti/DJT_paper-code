import sympy as sp
from sympy.physics.matrices import msigma

import partition
from DJT_matrix import *
import operators

N_c = operators.N_c
N_g = operators.N_g

q = 3/2
print("q =", q)
DJT = get_DJT(q)
DJT_dag = np.conj(DJT).T
U = operators.get_U(q = q)
L3 = operators.get_La(a=3, q=q, V=DJT, V_inv=DJT_dag)
N_alpha = partition.get_N_alpha(q)

q_max = q - 1/2 # only states with j up to (q - 1/2) satisfy the commutation relations
print("q_max =", q_max)

decimals = 15
print("rounding decimals", decimals)

print("The following loop tests the discrete eigenstates")
for j1 in [q_max - j_i/2 for j_i in range(0, int(2*q_max) + 1)]:
    print("-------------------------")
    print("  q =", sp.Rational(j1))
    deg_j1 = int(2*j1) + 1 # degeneracy of the left and right quantum numbers
    for m1 in [j1 - i_m for i_m in range(deg_j1)]:
        for mu1 in [j1 - i_mu for i_mu in range(deg_j1)]:
            print("(j, m, mu)=", sp.Rational(j1), sp.Rational(m1), sp.Rational(mu1))
            #
            v = get_eigenstate(j1, m1, mu1, q=q)
            for a in range(N_c):
                for b in range(N_c):
                    U_ab = operators.get_U_ab(U, a, b)
                    comm_ab = np.dot(L3, U_ab) - np.dot(U_ab, L3)
                    LHS = np.dot(comm_ab,v)
                    RHS = np.zeros(shape=(N_alpha, 1), dtype=complex)
                    for c in range(N_c):
                        tau3_ac = complex((msigma(3)[a,c]).evalf())/2.0
                        U_cb = operators.get_U_ab(U, c, b)
                        RHS = RHS + np.dot(tau3_ac*U_cb, v)
                    ####
                    msg = "(a,b)=({a}, {b}): |[L_3,U_ab]*v - \\tau^3_ac U_cb v|^2 = ".format(a=a, b=b)
                    diff = (LHS - RHS).round(decimals=decimals)
                    print(msg, get_norm2(diff))
                ####
            ####
            print("")
        ####
    ####
####
