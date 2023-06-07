import sympy as sp
##from sympy.physics.matrices import msigma

import partition
from DJT_matrix import *
import operators

N_c = operators.N_c
N_g = operators.N_g

q = 1
print("q =", q)
DJT = get_DJT(q)
DJT_dag = np.conj(DJT).T

U = operators.get_U(q = q)
U_dag = operators.get_Udag(U)

La = [operators.get_La(a=a, q=q, V=DJT, V_inv=DJT_dag) for a in [1, 2, 3]] 
Ra = [operators.get_Ra(a=a, q=q, V=DJT, V_inv=DJT_dag) for a in [1, 2, 3]] 
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
      v = get_DJT_column(DJT, j1, m1, mu1, q=q)
      for g in range(N_g):
        print("Generator: ", g)
        Lg = La[g]
        Rg = Ra[g]
        tau_g = operators.tau[g+1]
        for a in range(N_c):
          for b in range(N_c):
            U_ab = operators.get_U_ab(U, a, b)
            U_dag_ab = operators.get_U_ab(U_dag, a, b)
            #
            # [Lg, U] commutator
            comm1_ab = np.dot(Lg, U_ab) - np.dot(U_ab, Lg)
            comm2_ab = np.dot(Lg, U_dag_ab) - np.dot(U_dag_ab, Lg)
            LHS1 = np.dot(comm1_ab,v)
            LHS2 = np.dot(comm2_ab,v)
            RHS1 = np.zeros(shape=(N_alpha,), dtype=complex)
            RHS2 = np.zeros(shape=(N_alpha,), dtype=complex)
            for c in range(N_c):
              RHS1 = RHS1 + np.dot(-tau_g[a,c]*operators.get_U_ab(U, c, b), v)
              RHS2 = RHS2 + np.dot(operators.get_U_ab(U_dag, a, c)*tau_g[c,b], v)
            ##
            msg1 = "(a,b)=({a}, {b}): |[L_{g},U_ab]*v - \\sum_c (- \\tau^{g})_ac U_cb v|^2 = ".format(a=a, b=b, g=g+1) + 4*" "
            diff1 = (LHS1 - RHS1).round(decimals=decimals)
            ###print(msg1, get_norm2(diff1))
            msg2 = "(a,b)=({a}, {b}): |[L_{g},U^\dagger_ab]*v - \\sum_c U^\dagger_ac \\tau^{g}_cb v|^2 = ".format(a=a, b=b, g=g+1)
            diff2 = (LHS2 - RHS2).round(decimals=decimals)
            ###print(msg2, get_norm2(diff2))
            #
            # [Rg, U] commutator
            comm3_ab = np.dot(Rg, U_ab) - np.dot(U_ab, Rg)
            comm4_ab = np.dot(Rg, U_dag_ab) - np.dot(U_dag_ab, Rg)
            LHS3 = np.dot(comm3_ab,v)
            LHS4 = np.dot(comm4_ab,v)
            RHS3 = np.zeros(shape=(N_alpha,), dtype=complex)
            RHS4 = np.zeros(shape=(N_alpha,), dtype=complex)
            for c in range(N_c):
#              RHS3 = RHS3 + np.dot(operators.get_U_ab(U_dag, a, c)*tau_g[c,b], v)
              RHS3 = RHS3 + np.dot(-operators.get_U_ab(U,a,c)*tau_g[b,c], v)
              RHS4 = RHS4 + np.dot(tau_g[c,a]*operators.get_U_ab(U_dag, c, b), v)
            ##
            msg3 = "(a,b)=({a}, {b}): |[R_{g},U_ab]*v - \\sum_c (- \\tau^{g})_ac U_cb v|^2 = ".format(a=a, b=b, g=g+1) + 4*" "
            diff3 = (LHS3 - RHS3).round(decimals=decimals)
            print(msg3, get_norm2(diff3))
            msg4 = "(a,b)=({a}, {b}): |[R_{g},U^\dagger_ab]*v - \\sum_c U^\dagger_ac \\tau^{g}_cb v|^2 = ".format(a=a, b=b, g=g+1)
            diff4 = (LHS4 - RHS4).round(decimals=decimals)
            print(msg4, get_norm2(diff4))
          ####
        ####
        print("")
    ####
  ####
####
