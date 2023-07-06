# test the definition of operators on the spacetime lattice

import sys
sys.path.insert(0, '../')

import su2_DJT.lattice.operators as operators
from su2_DJT.lattice.grid import grid as grid

G = grid([3, 4, 2, 6])

# print(G.indices)

i = 15
print(i)
p = G.index_to_point(i)
print(p)
print(G.point_to_index(p))

print("ciao")
print(G.get_coordinates())
print("hello")
print(G.get_positions())


x, mu = G.index_to_x_mu(i)
i = G.x_mu_to_index(x=x, mu=mu)
x, mu = G.index_to_x_mu(i)
i = G.x_mu_to_index(x=x, mu=mu)

print(i)

##################

q = 1/2

G = grid([2, 3])
#U = operators.get_U(q=q, G=G)

#print(U)
