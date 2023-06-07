# import sympy as sp
# from DJT_matrix import *
# import operators

# q = 1
# DJT = get_DJT(q = q)

# q_max = q
# decimals = 15
# print("rounding decimals", decimals)

# print("The following loop tests the discrete eigenstates")
# for j1 in [q_max - j_i/2 for j_i in range(0, int(2*q_max) + 1)]:
#     print("-------------------------")
#     print("  q =", sp.Rational(j1))
#     deg_j1 = int(2*j1) + 1 # degeneracy of the left and right quantum numbers
#     for m1 in [j1 - i_m for i_m in range(deg_j1)]:
#         for mu1 in [j1 - i_mu for i_mu in range(deg_j1)]:
#             print("(j, m, mu)=", sp.Rational(j1), sp.Rational(m1), sp.Rational(mu1))
#             #
#             v = get_eigenstate(j1, m1, mu1, q=q)
#             print("|v|^2 =", get_norm2(v))
#             #
#             w1 = np.dot(L3, v) - m1*v
#             print("|L_3*v - m*v|^2 = ", get_norm2(w1).round(decimals=decimals))
#             #
#             w2 = np.dot(Lsquared, v) - j1*(j1+1)*v
#             print("|L^2*v - j*(j+1)*v|^2 = ", get_norm2(w2).round(decimals=decimals))
#             #
#             print("")
#         ####
#     ####
# ####
