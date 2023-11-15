# specturm of L^2: j(j+1)
# 2023-01/24 : I get the correct spectrum only for j even

from derivative import *
import matplotlib.pyplot as plt

n = 5

if n%2 == 0:
    print("Error: n should be odd. You gave:", n)
    quit()
####

g_m_max = int((n-1)/2) # psi \in [0, 4*pi]
print("g_m_max =", g_m_max)

print("Evaluating L^2")
Lsquared = get_Lsquared(n).evalf(chop=True)

print("Testing the spectrum")
n3 = int(n*n*n)
fun_0 = sp.zeros(n3, 1)

mu_dummy = 1 # dummy value, not relevant

for i_m_max in range(1, g_m_max+1):
    m_max = sp.Integer(i_m_max)/2
    print("===============")
    print("m_max =", m_max)
    print("===============")
    lambda_i = m_max*(m_max+1)
    for mz in range(-m_max, m_max+1):
        D = WignerD_vec(n, j=m_max, m=mz, mu=mu_dummy).evalf()
        # diff_eig = (Lsquared*D - lambda_i*D).expand(complex=True).evalf(12, chop=True)
        # print("Eigenfunction", mz, diff_eig==sp.zeros(n3, 1).evalf())
        eig = get_eig_matr_vec(Lsquared, D)
        print("Eigenvalues.", "Exact:", lambda_i, " - ", "Numerical:", eig)
    ####
####

