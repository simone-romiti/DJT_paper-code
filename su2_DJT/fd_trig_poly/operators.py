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

def get_Lplus(N, M):
    phi, theta = get_phi(N), get_theta(M)
    Dphi, Dtheta = get_Dphi(N), get_Dtheta(M)
    exp_iphi = apply_diag(np.exp, 1j*phi)
    cot_theta = apply_diag(np.cot, theta)
    res = tp(exp_iphi, Dtheta) + tp(exp_iphi*Dphi, 1j*cot_theta)
    return res
####

def get_Lminus(N, M):
    phi, theta = get_phi(N), get_theta(M)
    Dphi, Dtheta = get_Dphi(N), get_Dtheta(M)
    exp_miphi = apply_diag(np.exp, -1j*phi)
    cot_theta = apply_diag(np.cot, theta)
    res = tp(exp_miphi, -Dtheta) + tp(exp_miphi*Dphi, 1j*cot_theta)
    return res.evalf(chop=True)
####

def get_Lx(N, M):
    Lplus, Lminus = get_Lplus(N, M), get_Lminus(N, M)
    Lx = (Lplus+Lminus)/2
    return Lx.evalf(chop=True)
####

def get_Ly(N, M):
    Lplus, Lminus = get_Lplus(N, M), get_Lminus(N, M)
    Ly = (Lplus-Lminus)/(2*1j)
    return Ly.evalf(chop=True)
####


def get_Lz(N):
    Lz = -1j*get_Dphi(N)
    return Lz.evalf()
####

def get_L2(N, M):
    theta = get_theta(M)
    Dphi = get_Dphi(N)
    Dphi2 = Dphi*Dphi
    Dtheta = get_Dtheta(M)
    Dtheta2 = Dtheta*Dtheta
    cot_theta = apply_diag(lambda z: 1/np.tan(z), theta)
    sin_m2_theta = apply_diag(sin_m2, theta)
    #
    IdN = np.matrix(np.eye(N))
    T1 = np.matrix(Dtheta2 + cot_theta*Dtheta)
    res = np.kron(IdN, T1) + np.kron(Dphi2, sin_m2_theta)
    res = -res
    return np.matrix(res)
####

def get_fourier_mode(N, s):
    Ui = np.matrix(np.zeros(N,1))
    for i in range(N):
        phi_i = i*(2*np.pi/N)
        Ui[i, 0] = sp.exp(s*1j*phi_i)
    ####
    return Ui
####

