# operators.py
# Reference: https://arxiv.org/pdf/quant-ph/0008120.pdf

import numpy as np
import sympy as sp
from sympy.physics.quantum import TensorProduct as tp
import matplotlib.pyplot as plt

def get_t(x, xl):
    N = len(xl)
    t = 1
    for l in range(N):
        t *= sp.sin((x-xl[l])/2)
    ####
    return t
####

def get_tprime(x, xl):
    t = get_t(x, xl)
    return sp.diff(t, x)
#####


def get_T(N, xl):
    x = sp.symbols("x")
    tprime = get_tprime(x, xl)
    T = np.matrix(np.zeros(shape=(N,N)))
    for i in range(N):
        T[i,i] = tprime.subs(x, xl[i]).evalf()
    ####
    return T
####

def get_T_inv(N, xl):
    T = get_T(N, xl)
    T_inv =np.zeros(shape=(N,N))
    for i in range(N):
        T_inv[i,i] = 1/T[i,i]
    ####
    return T_inv
####

def get_Dtilde(N, x):
    Dtilde = np.matrix(np.zeros(shape=(N,N)))
    for i in range(N):
        for j in range(N):
            if i==j:
                s = 0
                for l in range(N):
                    if i==l:
                        continue
                    else:
                        s += (1/2)*1/np.tan((x[i]-x[l])/2)
                    ####
                ####
                Dtilde[i,j] = s
            else:
                Dtilde[i,j] = (1/2)*1/np.sin((x[i]-x[j])/2)
            ####
        ####
    ####
    return Dtilde
####

def get_D_nodes(N, xj):
    Dtilde = get_Dtilde(N, xj)
    T = get_T(N, xj)
    T_inv = get_T_inv(N, xj)
    return T*Dtilde*T_inv
####

def get_nodes(N, xmin, xmax):
    delta = (xmax-xmin)/N
    return [xmin + delta*j for j in range(1, N+1)]
####

def get_D(N, xmin, xmax):
    xj = get_nodes(N, xmin, xmax)
    return get_D_nodes(N, xj)
####

def get_angle(N, a1, a2):
    alpha = np.matrix(np.zeros(shape=(N,N)))
    alpha_vals = get_nodes(N, a1, a2)
    for i in range(N):
        alpha[i,i] = alpha_vals[i]
    ####
    return alpha
####

def get_phi(N):
    return get_angle(N, -np.pi, np.pi)

def get_Dphi(N):
    return get_D(N, -np.pi, np.pi)
####

def get_Dpsi(N):
    return get_D(N, -np.pi, np.pi)
####

def get_theta(M):
    return get_angle(M, 0, np.pi - 1/M)
####

def get_Dtheta(M):
    return get_D(M, 0, np.pi - 1/M)
####

def apply_diag(fun, x):
    """ apply function to diagonal matrix F_{ij} = fun(x_{ij}) """
    n = x.shape[0]
    fx = np.matrix(np.zeros(shape=(n,n)))
    for i in range(n):
        fx[i,i] = fun(x[i, i])
    ####
    return fx
####

def sin_m2(x):
    return np.sin(x)**(-2)
####

class momenta:
    def __init__(self, theta, phi, psi):
        self.theta = theta
        self.phi = phi
        self.psi = psi
        self.N_theta = len(list(self.theta))
        self.N_phi = len(list(self.phi))
        self.N_psi = len(list(self.psi))
        self.n = self.N_theta * self.N_phi * self.N_psi

        ## partial derivatives
        self.D_theta = get_Dtheta(self.N_theta)
        self.D_phi = get_Dphi(self.N_phi)
        self.D_psi = get_Dpsi(self.N_psi)

        Id_theta = np.matrix(np.eye(self.N_theta))
        Id_phi = np.matrix(np.eye(self.N_phi))
        Id_psi = np.matrix(np.eye(self.N_psi))
        self.kron_D_theta = np.kron(np.kron(self.D_theta, Id_phi), Id_psi)
        self.kron_D_phi = np.kron(np.kron(Id_theta, self.D_phi), Id_psi)
        self.kron_D_psi = np.kron(np.kron(Id_theta, Id_phi), self.D_psi)

        ## initializing some useful matrices
        self.cot_theta = np.kron(np.kron(apply_diag(lambda t: 1/np.tan(t), get_theta(self.N_theta)), Id_phi), Id_psi)
        self.cos_phi = np.kron(np.kron(Id_theta, apply_diag(np.cos, get_phi(self.N_phi))), Id_psi)
        self.sin_phi = np.kron(np.kron(Id_theta, apply_diag(np.sin, get_phi(self.N_phi))), Id_phi)
        self.cos_phi = np.kron(np.kron(Id_theta, apply_diag(np.cos, get_phi(self.N_phi))), Id_phi)
        self.inv_sin_theta = np.kron(np.kron(apply_diag(lambda t: 1/np.sin(t), get_theta(self.N_theta)), Id_phi), Id_psi)
    ####
    def Lsquared(self):
        A1 = -(self.cot_theta)*self.kron_D_theta
        A2 = -(self.kron_D_theta**2)
        A3 = -(self.cot_theta**2)*(self.kron_D_phi**2) - (self.kron_D_phi**2)
        A4 = +2*(self.cot_theta)*(self.inv_sin_theta)*(self.kron_D_phi*self.kron_D_psi)
        A5 = - (self.inv_sin_theta**2)*(self.kron_D_psi**2)
        return A1+A2+A3+A4+A5
    ####
    def L1(self):
        a = - self.cot_theta * self.cos_phi * self.D_phi
        b = - self.sin_phi * self.D_theta
        c = + self.cos_phi * self.inv_sin_theta * self.D_psi
        res = -1j * (a+b+c)
        return res
    ####
    def L2(self):
        a = - self.cot_theta * self.sin_phi * self.D_phi
        b = + self.cos_phi * self.D_theta
        c = + self.sin_phi * self.inv_sin_theta * self.D_psi
        res = -1j * (a+b+c)
        return res
    ####
    def L3(self):
        return -1j * self.D_phi
    ####
    def Lplus(self):
        return self.L1() + 1j * self.L2()
    ####
    def Lminus(self):
        return self.Lplus().getH()
    ####
####

