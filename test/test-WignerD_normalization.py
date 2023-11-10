import sympy as sp
import numpy as np
from su2_DJT.DJT.DJT_matrix import get_WignerD

onehalf = (1/2)
theta, phi, psi = 0.1, -0.2, -0.3

x1 = get_WignerD(onehalf, -onehalf, -onehalf, phi, theta, psi, use_sympy=False)
x2 = np.cos(theta/2)*np.exp(-1j*(phi+psi)/2)
print(x1, x2, x1 - x2)

y1 = - get_WignerD(onehalf, -onehalf, +onehalf, phi, theta, psi, use_sympy=False)
y2 = - np.sin(theta/2)*np.exp(-1j*(phi-psi)/2)
print(y1, y2, y1 - y2)

z1 = - get_WignerD(onehalf, +onehalf, -onehalf, phi, theta, psi, use_sympy=False)
z2 = + np.sin(theta/2)*np.exp(+1j*(phi-psi)/2)
print(z1, z2, z1 - z2)


w1 = get_WignerD(onehalf, +onehalf, +onehalf, phi, theta, psi, use_sympy=False)
w2 = + np.cos(theta/2)*np.exp(+1j*(phi+psi)/2)
print(w1, w2, w1 - w2)
