# testing hermiticity on eigenfunctions

from derivative import *
import matplotlib.pyplot as plt

n = 7

if n%2 == 0:
    print("Error: n should be odd. You gave:", n)
#    quit()
####

g_m_max = int((n-1)/2) # psi \in [0, 4*pi]
print("g_m_max =", g_m_max)

print("Evaluating the matrices")
L3 = get_L3(n).evalf(chop=True).H
Lplus = get_Lplus(n).evalf(chop=True)
Lminus = Lplus.H # get_Lminus(n).evalf(chop=True)

print("Testing the commutation relations")
n3 = int(n*n*n)
fun_0 = sp.zeros(n3, 1)

mu_dummy = 0 # dummy value, not relevant

for jp in range(0, 2*g_m_max+1):
    m_max = sp.Integer(jp)/2
    for nz in range(0, jp+1):
        mz = sp.Rational(-m_max + nz)
        print("j=",m_max, "mz=", mz)
        D = WignerD_vec(n, j=m_max, m=mz, mu=mu_dummy).evalf()
        diff_eig = (get_eig_matr_vec(L3, D) - mz).expand(complex=True).evalf(12, chop=True)
        print("Eigenfunction", mz, diff_eig==0)
        Lplus_fact = -sp.sqrt(m_max*(m_max+1) - mz*(mz+1)).evalf(chop=True)
        Dplus = WignerD_vec(n, j=m_max, m=mz+1, mu=mu_dummy).evalf(chop=True)
        diff_plus = (Dplus.H*Lplus*D - Lplus_fact*Dplus.H*Dplus).evalf(12,chop=True)
        print("Ladder Lplus:", (Dplus.H*Lplus*D).evalf(chop=True)[0,0], (Lplus_fact*Dplus.H*Dplus).evalf(12, chop=True)[0,0])
        Lminus_fact = sp.sqrt(m_max*(m_max+1) - mz*(mz-1)).evalf()
        Dminus = WignerD_vec(n, j=m_max, m=mz-1, mu=mu_dummy).evalf()
        diff_minus = (Dminus.H*Lminus*D - Lminus_fact*Dminus.H*Dminus).evalf(12,chop=True)[0,0]
        print("Ladder Lminus:", diff_minus==0)
    ####
####

