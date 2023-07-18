# testing the size of the matrix elements of the momenta around U=1

import sys
sys.path.insert(0, '../')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle

import su2_DJT.operators.operators as operators
from   su2_DJT.DJT.DJT_matrix import *
import su2_DJT.S3_sphere.partition as partition
import su2_DJT.S3_sphere.indices as indices

for q in [1/2, 2/2, 3/2, 4/2, 5/2, 6/2]:
    DJT = get_DJT(q)
    DJT_dag = get_DJT_dag(DJT = DJT)

    Lsquared = operators.get_Lsquared(q=q, V=DJT, V_inv=DJT_dag)

    N_alpha = Lsquared.shape[0]

    x, y = [], []
    for i in range(N_alpha):
        alpha_i = indices.S3_point_to_angles_value(i, q=q)
        y_i = partition.get_yi(alpha_i)
        for j in range(N_alpha):
            alpha_j = indices.S3_point_to_angles_value(j, q=q)
            y_j = partition.get_yi(alpha_j)
            d_ij = partition.get_distance_S3(y_i, y_j)
            x.append(d_ij)
            y.append(np.abs(Lsquared[i,j]))
        ####
    ####

    plt.plot(x, y, linestyle="None", marker=".")
    plt.savefig("sphere/distance_VS_value-q{q}.pdf".format(q=q))
    plt.cla()

    plt.hist(x)
    plt.savefig("sphere/distances-q{q}.pdf".format(q=q))
    plt.cla()

    plt.hist(y)
    plt.savefig("sphere/values-q{q}.pdf".format(q=q))
    plt.cla()
####


