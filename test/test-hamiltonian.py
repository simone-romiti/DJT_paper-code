# test the definition of the Hamiltonian

import sys
sys.path.insert(0, '../')

import su2_DJT.lattice.hamiltonian as hamiltonian
from su2_DJT.lattice.grid import grid as grid
from su2_DJT.DJT.DJT_matrix import *

G = grid([3, 3, 2], 3)

q = 1/2
DJT = get_DJT(q)
DJT_dag = get_DJT_dag(DJT = DJT)

g = 0.1

#H_mag = hamiltonian.get_H_mag(g, G, q)

H_el = hamiltonian.get_H_el(g, G, q, V=DJT, V_inv=DJT_dag)

