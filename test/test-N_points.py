# this script shows how N_alpha - N_q is always even

import sympy as sp
import numpy as np


import sys
sys.path.insert(0, '../')

import su2_DJT.S3_sphere.partition as partition

q_max = 30

for q in [x/2  for x in range(int(2*q_max))]:
    #print("q =", q)
    N_q = partition.get_N_q(q)
    N_alpha = partition.get_N_alpha(q)
    dN = N_alpha-N_q
    if (dN % 2) != 0:
        print("not even")
    ####    
####
