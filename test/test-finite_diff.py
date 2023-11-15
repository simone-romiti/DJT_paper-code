
import numpy as np
decimals = 14

import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '../')

from su2_DJT.fd_trig_poly.operators import *
from su2_DJT.S3_sphere.partition import get_N_q

"""
This is the setup that works:

q_max integer:
N1 = int(2*q_max + 1)
N2 = int(4*q_max + 1)
N3 = int(4*q_max + 1)

q_max half-integer:
N1 = int(2*q_max + 1)
N2 = int(4*q_max + 2)
N3 = int(4*q_max + 2)
"""

q_max = 3/2
N1 = int(2*q_max + 1)
N2 = int(4*q_max + 2)
N3 = int(4*q_max + 2)

if q_max.is_integer():
    N_qmax = int(sum([(2*j + 1)**2 for j in range(q_max+1)]))
else:
    N_qmax = int(sum([(j + 1)**2 for j in range(-int(2*q_max), int(2*q_max)+1, 2)]))
####

theta = [k*np.pi/(N1+1) for k in range(1, N1+1)]
phi = [k*4*np.pi/(N2+1) for k in range(1, N2+1)]
psi = [k*4*np.pi/(N3+1) for k in range(1, N2+1)]

mom = momenta(theta, phi, psi)
L2 = mom.Lsquared()

eigs = np.real(np.linalg.eigvals(L2))
eigs = [x for x in sorted(eigs)]
eigs = eigs[0:N_qmax]

plt.plot(eigs, linestyle="None", marker="*", label="continuum")
plt.show()
