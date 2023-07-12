# distribution of the matrix elements on the sphere S_3


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys
sys.path.insert(0, '/home/simone/Documents/su2_discrete/su2_discrete-code/')

import su2_DJT.S3_sphere.partition as partition
import su2_DJT.S3_sphere.indices as indices

from su2_DJT.operators.operators import *
from su2_DJT.DJT.DJT_matrix import get_norm, get_DJT, get_DJT_dag

# distribution of elements with respect to 1
def distribution_wrt_1(M, q):
    n = M.shape[0]
    X, Y = [], []
    alpha_0 = [0.0, 0.0, 0.0]
    y_0 = np.matrix(np.array(partition.get_yi(alpha_0)))
    idx = 0
    for i in range(n):
#        alpha_i = indices.S3_point_to_angles_value(i, q)
#        y_i = np.matrix(np.array(partition.get_yi(alpha_i)))
#        y_i = np.reshape(y_i, newshape=(4,1))
        for j in range(n):
#            alpha_j = indices.S3_point_to_angles_value(j, q)
#            y_j = np.matrix(np.array(partition.get_yi(alpha_j)))
#            y_j = np.reshape(y_j, newshape=(4,1))
            #
            alpha_idx = indices.S3_point_to_angles_value(idx, q)
            y_idx = np.matrix(np.array(partition.get_yi(alpha_idx)))
            y_idx = np.reshape(y_idx, newshape=(4,1))
            dy_norm = get_norm(y_idx-y_0)
            X.append(dy_norm)
            Y.append(np.absolute(M[i,j]))

            idx += 1
        ####
    ####
    return [X, Y]
####

q = 1
DJT = get_DJT(q=q)
DJT_dag = get_DJT_dag(DJT=DJT)
Lsquared = get_Lsquared(q=q, V=DJT, V_inv=DJT_dag)
#U_00 = get_U(q=q)[0][0]

D = distribution_wrt_1(Lsquared, q)

x, y = D
x, y = zip(*sorted(zip(x, y)))

#counts, bins = np.histogram(y, density=True)
#plt.stairs(counts, bins)

fig, ax = plt.subplots()

n_bins = partition.get_N_phi(q=q) * partition.get_N_psi(q=q)
n_bins = int(np.sqrt(n_bins))
#counts, xedges, yedges, im = ax.hist2d(x, y, bins=n_bins, norm=LogNorm())
#fig.colorbar(im, ax=ax)

#plt.hist2d(x, y, bins=int(np.sqrt(len(x))))
plt.xlabel("d")
plt.ylabel("|M_ij|")
# plt.legend()

plt.plot(x, y)

plt.savefig("plot.pdf")
#plt.show()


