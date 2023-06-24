import numpy as np
import sympy as sp
from sympy.physics.quantum.spin import WignerD

import sys
sys.path.insert(0, '../')

import su2_DJT.S3_sphere.partition as partition
import su2_DJT.S3_sphere.indices as indices
import su2_DJT.operators.operators as operators


q = 3/2
print("q=", q)

decimals = 15
print("rounding decimals", decimals)

Uc = operators.get_U(q)
Uc_dag = operators.get_Udag(Uc)

N_c = operators.N_c

N_alpha = partition.get_N_alpha(q=q)
N_theta = partition.get_N_theta(q)
N_phi = partition.get_N_phi(q)
N_psi = partition.get_N_psi(q)
N_alpha = partition.get_N_alpha(q)
theta = partition.get_theta(N_theta)
phi = partition.get_phi(N_phi)
psi = partition.get_psi(N_psi)

onehalf = sp.Rational(1/2)
for i in range(N_alpha):
  i_theta, i_psi, i_phi = indices.S3_point_to_angles_index(i, q)
  W = [[0,0],[0,0]]
  W[0][0] = complex(+WignerD(onehalf, +onehalf , +onehalf, phi[i_phi], theta[i_theta], psi[i_psi]).doit().evalf())
  W[0][1] = complex(-WignerD(onehalf, -onehalf , +onehalf, phi[i_phi], theta[i_theta], psi[i_psi]).doit().evalf())
  W[1][0] = complex(-WignerD(onehalf, -onehalf , +onehalf, phi[i_phi], theta[i_theta], psi[i_psi]).doit().evalf())
  W[1][1] = complex(+WignerD(onehalf, -onehalf , -onehalf, phi[i_phi], theta[i_theta], psi[i_psi]).doit().evalf())
  # w00 = sp.cos(theta[i_theta]/2)*sp.exp(-sp.I*(psi[i_psi]+phi[i_phi])/2).evalf()
  w01 = complex(-sp.sin(theta[i_theta]/2)*sp.exp(sp.I*(-psi[i_psi]+phi[i_phi])/2).evalf(chop=True))
#  print(WignerD(onehalf, +onehalf , -onehalf, phi[i_phi], theta[i_theta], psi[i_psi]).doit().evalf() + np.conj(WignerD(onehalf, -onehalf , +onehalf, phi[i_phi], theta[i_theta], psi[i_psi]).doit().evalf()))
  print(np.complex64(W[0][1] - w01).round(decimals=decimals))
  for a in range(N_c):
    for b in range(N_c):
      U = Uc[a][b]
      U_dag = Uc_dag[a][b]
      # print(a, b, i, U[i,i], U_dag[i,i],  W[a][b]) 
#    print(a, b, i, U_dag[i,i],  W[a][b]) 
    ##
  ##
##

