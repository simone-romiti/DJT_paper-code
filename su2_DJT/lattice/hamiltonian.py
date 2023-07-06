# electric and magnetic Hamiltonian contributions

import numpy as np

import su2_DJT.S3_sphere.partition as partition
import su2_DJT.lattice.operators as operators

def get_H_mag(g, G, q):
    points = G.get_positions()
    d = G.get_d()
    N_p = G.get_N() # total number of pairs (x, mu)
    N_alpha = partition.get_N_alpha(q=q)
    N = N_p*N_alpha
    H_mag = np.matrix(np.zeros(shape=(N, N), dtype=complex))
    for x in points:
        for mu in range(d):
            for nu in range(mu+1, d):
                H_mag += operators.get_plaquette(G=G, x=x, mu=mu, nu=nu, q=q)
            ####
        ####
    ####
    H_mag = H_mag + H_mag.getH()
    H_mag = (1/(4*g*g))*H_mag
    return H_mag
####

def get_H_el(g, G, q, V, V_inv):
    points = G.get_positions()
    d = G.get_d()
    N_p = G.get_N() # total number of pairs (x, mu)
    N_alpha = partition.get_N_alpha(q=q)
    N = N_p*N_alpha
    H_el = np.matrix(np.zeros(shape=(N, N), dtype=complex))
    for x in points:
        for mu in range(d):
            H_el += operators.get_Lsquared(G=G, x=x, mu=mu, q=q, V=V, V_inv=V_inv)
        ####
    ####
    H_el = (g*g/2)*H_el
    return H_el
####

