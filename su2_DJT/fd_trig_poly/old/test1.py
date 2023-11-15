

from derivative import *

n = 5
m_max = int((n-1)/2)

T = 2*sp.pi
al = get_nodes_vals(n, 0, T)
Lz = -sp.I*get_Dangle(n, T, al)

for mz in range(-m_max, m_max+1):
    f1 = fun1_to_vec(sp.exp, sp.I*mz*al)
    diff = (Lz*f1 - mz*f1).evalf(chop=True)
    print(mz, diff)
####

