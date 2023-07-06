# lattice grid of points in d dimensions
from functools import reduce
import operator
import numpy as np

class grid:
    def __init__(self, Ni):
        self.Ni = Ni # dimension sizes
        self.d = len([x for x in Ni if x>1]) # number of dimensions
        self.Npts = reduce(operator.mul, self.Ni) # total number of points
        idx_array = np.arange(self.Npts)
        idx_array = np.reshape(idx_array, newshape=Ni)
        self.index_enum = [(idx, p) for idx, p in np.ndenumerate(idx_array)]
        self.indices = [idx for idx, _ in self.index_enum]
        self.coordinates = [p for _, p in self.index_enum]
    ####
    def get_Npts(self):
        return self.Npts
    ####
    def get_sizes(self):
        return self.Ni
    ####
    def get_index_enum(self):
        return self.index_enum
    ####
    def index_to_point(self, i):
        return self.coordinates[i]
    ####
    def point_to_index(self, p):
        c = np.array(self.coordinates)
        idx = np.where(np.equal(c, p))
        return idx[0][0]
    ####        
####

G = grid([3, 4, 2, 6])

# print(G.indices)

i = 15
print(i)
p = G.index_to_point(i)
print(p)
print(G.point_to_index(p))



