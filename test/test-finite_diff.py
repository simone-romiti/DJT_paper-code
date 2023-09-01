
import numpy as np
decimals = 14

import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '../')

from su2_DJT.fd_trig_poly.operators import *
#from su2_DJT.fd_trig_poly.derivative import *


n = 3 # 2*lmax +1
lmax = (n-1)/2
print("lmax:" , lmax)

L2 = np.matrix(get_L2(n, 2*n+1), dtype=complex)

eigs = np.real(np.linalg.eigvals(L2))
eigs = [x for x in sorted(eigs)]

plt.plot(eigs, linestyle="None", marker="*", label="continuum")
plt.show()

# #import su2_DJT.fd_trig_poly.derivative as fd

# for N in [14, 15]:
#     tau = 4*np.pi
#     alpha = np.array([k*tau/(N+1) for k in range(1, N+1)])
#     D = -1j * D_alpha(tau=tau, x=alpha).get_matrix() ## get_D(N)
#     eigs = np.linalg.eigvals(D)
#     eigs = np.round(eigs, decimals=decimals)
#     print(eigs)
#     m  = 1/2
#     f1 = get_function(lambda t: np.exp(1j * m * t), alpha)
#     f2 = get_function(lambda t: np.exp(1j * m * t), alpha)
#     print(np.array(D*f1)/np.array(m*f2))
# ####

# for N in [15]:
#     tau = 4*np.pi
#     alpha = np.array([k*tau/(N+1) for k in range(1, N+1)])
#     D = -1j * (2*np.pi/tau) *D_alpha(tau=tau, x=alpha).get_matrix() ## get_D(N)
#     eigs = np.linalg.eigvals(D)
#     eigs = np.round(eigs, decimals=decimals)
#     print(eigs)
#     m  = 1/2
#     f1 = get_function(lambda t: np.exp(1j*m*t), alpha)
#     f2 = get_function(lambda t: np.exp(1j*m*t), alpha)
#     print(np.array(D*f1)/np.array(m*f2))
# ####