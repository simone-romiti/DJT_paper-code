
import sympy as sp
from DJT_matrix import *

q = 2

print("The following output should be a list of all True")

for j1 in [q - j_i/2 for j_i in range(0, int(2*q) + 1)]:
    deg_j1 = int(2*j1) + 1 # degeneracy of the left and right quantum numbers
    for m1 in [j1 - i_m for i_m in range(deg_j1)]:
        for mu1 in [j1 - i_mu for i_mu in range(deg_j1)]:
            ##print(j1, m1, mu1)
            i1 = su2_irrep_to_index(j1, m1, mu1, q)
            j2, m2, mu2 = [float(x) for x in su2_index_to_irrep(i1, q)]
            print([j1, m1, mu1] == [j2, m2, mu2])
        ####
    ####
####
