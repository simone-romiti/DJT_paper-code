
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.size": 18
})

import sys
sys.path.insert(0, '../')

from su2_DJT.S3_sphere.partition import get_N_q
from su2_DJT.DJT.DJT_matrix import get_DJT, get_DJT_dag, get_P_garbage_space
import su2_DJT.operators.operators as operators

q = 1.5
N_q = get_N_q(q=q)
DJT = get_DJT(q = q)
DJT_dag = get_DJT_dag(DJT)
Lsquared = operators.get_Lsquared(q=q, V=DJT, V_inv=DJT_dag)

q_proj = q + 10
kappa = q_proj*(q_proj+1)
P = get_P_garbage_space(DJT, q)
Lsquared = Lsquared + kappa*P

# Notes: 
# - the eigenvalues are real by construction
# - N_\alpha - N_q eigenstates are projected above q*(q+1)

Lsqr_eigs, _ = np.linalg.eig(Lsquared)
Lsqr_eigs = [l_i for l_i in Lsqr_eigs if l_i < (q+1)*(q+2)]
Lsqr_eigs = sorted(np.real(Lsqr_eigs))

# Plot the eigenvalues
plt.scatter(
    [i for i in range(len(Lsqr_eigs))], Lsqr_eigs, 
    linestyle="None", marker="x", s=20)
plt.xlabel('$i$')
plt.ylabel('$\lambda_i$')
#plt.title('Eigenvalues of $\sum_a L_a L_a$')
plt.grid(True, linestyle="dashed")
plt.tight_layout()
plt.savefig("./Lsquared.pdf")

plt.close()


markers = ["1", "2", "3"]
for a in range(3):
    La = operators.get_La(a=a+1, q=q, V=DJT, V_inv=DJT_dag)
    La = La + kappa*P
    La_eigs, _ = np.linalg.eig(La)
    La_eigs += (a-1)*2.5e-2
    La_eigs = [l_i for l_i in La_eigs if l_i < (q+1)*(q+2)]
    La_eigs = sorted(np.real(La_eigs))

    # Plot the eigenvalues
    plt.scatter(
        [i for i in range(len(La_eigs))], La_eigs, 
        linestyle="None", marker=markers[a], s=20,
        label = "$L_"+str(a+1)+"$")
####
plt.xlabel('$i$')
plt.ylabel('$\lambda_i$')
#plt.title('Eigenvalues of $\sum_a L_a L_a$')
plt.grid(True, linestyle="dashed")
plt.legend()
plt.tight_layout()
plt.savefig("./La.pdf")

