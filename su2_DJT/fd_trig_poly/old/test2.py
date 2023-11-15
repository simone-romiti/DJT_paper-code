# testing action of L3, Lplus and Lminus on eigenfunctions

from derivative import *
import matplotlib.pyplot as plt

n = 9

if n%2 == 0:
    print("Error: n should be odd. You gave:", n)
    quit()
####

g_m_max = int((n-1)/2) # psi \in [0, 4*pi]
print("g_m_max =", g_m_max)

print("Evaluating the matrices")
L3 = get_L3(n).evalf(chop=True)
Lplus = get_Lplus(n).evalf(chop=True)
Lminus = get_Lminus(n).evalf(chop=True)

print("Testing the commutation relations")
n3 = int(n*n*n)
fun_0 = sp.zeros(n3, 1)

mu_dummy = 0 # dummy value, not relevant

for m_max in range(g_m_max, g_m_max+1):
    print("===============")
    print("m_max =", m_max)
    print("===============")
    for mz in range(-m_max, m_max+1):
        D = WignerD_vec(n, j=m_max, m=mz, mu=mu_dummy).evalf()
        diff_eig = (L3*D - mz*D).expand(complex=True).evalf(12, chop=True)
        print("Eigenfunction", mz, diff_eig==sp.zeros(n3, 1).evalf())
        diff_plus = (L3*Lplus*D - (mz+1)*Lplus*D).evalf(12,chop=True)
        print("Ladder Lplus:", diff_plus==fun_0)
        diff_minus = (L3*Lminus*D - (mz-1)*Lminus*D).evalf(12,chop=True)
        print("Ladder Lminus:", diff_minus==fun_0)
    ####
####

