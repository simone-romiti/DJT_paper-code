# lattice grid of points in d dimensions
from functools import reduce
import operator
import numpy as np

def get_index_enum(N, sizes):
    idx_array = np.arange(N)
    idx_array = np.reshape(idx_array, newshape = sizes)
    index_enum = [(idx, p) for idx, p in np.ndenumerate(idx_array)]
    return index_enum
####

def  get_indices_and_coordinates(index_enum):
    indices = [idx for _, idx in index_enum]
    coordinates = [p for p, _ in index_enum]
    return (indices, coordinates)
####

class grid:
    def __init__(self, Ni, d):
        self.Ni = Ni # dimension sizes
        self.d = d
        if(len(Ni) != d):
            print("sizes :", Ni)
            print("d : ", d)
            raise ValueError("Invalid lattice grid")
        ####
        self.N_pos = reduce(operator.mul, self.Ni) # total number of points
        self.N = self.N_pos*self.d # one operator for each position AND direction
        self.index_enum = get_index_enum(self.N, sizes = list(Ni)+[self.d])
        self.indices, self.coordinates = get_indices_and_coordinates(self.index_enum)
    ####
    def get_N_pos(self):
        return self.N_pos
    ####
    def get_N(self):
        return self.N
    ####
    def get_sizes(self):
        return self.Ni
    ####
    def get_d(self):
        return self.d
    ####
    def get_index_enum(self):
        return self.index_enum
    ####
    def get_coordinates(self):
        return self.coordinates
    ####
    def get_positions(self):
        P = get_indices_and_coordinates(get_index_enum(self.N_pos, sizes = self.Ni))[1]
        return P
    ####
    def index_to_point(self, i):
        return self.coordinates[i]
    ####
    def point_to_index(self, p):
        c = np.array(self.coordinates)
        idx = np.where(np.all(c == p, axis=1))
        return idx[0][0]
    ####
    def index_to_x_mu(self, i):
        c = self.index_to_point(i=i)
        return (list(c[:-1]), c[-1]) # (x, mu) 
    ####
    def x_mu_to_index(self, x, mu):
        p = list(x) + [mu]
        return self.point_to_index(p)
    ####        
####

