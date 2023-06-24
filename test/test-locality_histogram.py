# testing the locality of the operator whan the truncation approached infinity


import sys
sys.path.insert(0, '../')

import su2_DJT.operators.operators as operators
from   su2_DJT.DJT.DJT_matrix import *
#import S3_sphere.partition as partition

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle


q_vals = [k/2 for k in range(1, 6)]
q_max = max(q_vals)
N_g = operators.N_g # number of generators in the algebra

abs_La_dict = {1: [], 2: [], 3: []}

print("Building the matrices for the L_a")
for q in q_vals:
    print("q =", q)
    DJT = get_DJT(q)
    DJT_dag = get_DJT_dag(DJT = DJT)
    #
    for a in range(1, N_g+1):
        La = operators.get_La(a=a, q=q, V=DJT, V_inv=DJT_dag)
        abs_values = np.abs(np.array(La))  # Compute absolute values of elements
        flat_values = abs_values.flatten()
        N = flat_values.shape[0]
        plt.hist(x=flat_values, bins=int(np.sqrt(N)), density=True)
#        plt.legend()
        plt.savefig("./histo/q{q}_a{a}.pdf".format(q=q, a=a) )
        plt.cla()
    ####
####

