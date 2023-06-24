# operators defined on a spatial grid of points

import su2_DJT.operators as operators

import numpy as np
import functools

# lattice grid in "d" dimensions
class lattice:
    def __init__(self, Ni):
        self.d = len(Ni) # number of dimensions
        self.Li = Ni # list of number of points along each direction
        self.V = functools.reduce(lambda x, y: x * y, Ni)  # total number of points
        self.N = self.V * self.d # number of gauge links
    ####
    def get_N(self):
        return self.N
    ####
####



## U(i) = U(x, mu) for a given pair of (x, mu) in the lattice Lat
def get_U(Lat: lattice, i: int):
    N = Lat.get_N()
    U = np.matrix(np.zeros(shape=(N,N)))
    for i in range(N):
        U[]


    return U
####


