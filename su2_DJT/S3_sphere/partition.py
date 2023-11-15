""" Polynomial partitioning of the sphere """

import numpy as np
import sympy as sp

def is_half_integer(j):
    b1 = (float(j)).is_integer()
    b2 = (float(2.0*j)).is_integer()
    return (not b1) and b2
####

def check_truncation(q):
    b = (float(2.0*q)).is_integer()
    if not b:
        raise ValueError("Illegal truncation: 2*{q} is an integer".format(q=q))
    ####
####

def get_N_q(q):
    """ \sum_{j=0}^{q} \sum_{m, \mu = -j}^{j}  (2j +  1)^2 """
    check_truncation(q)
    N_q = (4*q + 3)*(2*q + 2)*(2*q + 1) # is always a multiple of 6
    N_q = int(N_q/6)
    return N_q
####

def get_N_theta(q):
    delta = 1
    if not (q+0.0).is_integer():
        delta = 1/2
    ####
    return int(q+delta)
####

def get_N_phi(q):
    return int(4*q) + 1
####

def get_N_psi(q):
    return int(4*q) + 1
####

def get_N_alpha(q):
    check_truncation(q)
    N_theta = get_N_theta(q)
    N_phi = get_N_phi(q)
    N_psi = get_N_psi(q)
    N_alpha = N_theta*N_phi*N_psi
    return N_alpha
####


def get_asympt_dtheta(q):
    """ infinitesimal increment in theta when q \to \infty """
    return np.pi/get_N_theta(q=q)
####

def get_dphi(q):
    """ increment in phi """
    return 4*np.pi/get_N_phi(q=q)
####

def get_dpsi(q):
    """ increment in psi """
    return 4*np.pi/get_N_psi(q=q)
####

def get_asympt_d3alpha(q):
    """ measure of the infinitesimal increment for q \to \infty """
    dtheta = get_asympt_dtheta(q)
    dphi = get_dphi(q)
    dpsi = get_dpsi(q)
    return dtheta*dphi*dpsi
####

def get_Legende_roots_x(n):
    """ roots of the n-th Legendre polynomial: P_n(x) """
    # Calculate the coefficients of the n-th Legendre polynomial
    coeffs = np.zeros(n + 1)
    coeffs[n] = 1
    legendre_poly = np.polynomial.legendre.Legendre(coeffs)

    # Find the roots of the Legendre polynomial
    roots = legendre_poly.roots()

    return roots
####

def get_Legendre_roots_theta(n):
    """ roots of the n-th Legendre polynomial P_n(x = \cos(theta)) """
    roots = get_Legende_roots_x(n)
    return np.arccos(roots)
####

def get_theta(N_theta):
    theta = get_Legendre_roots_theta(N_theta)
    return np.array(theta)
####

def get_phi(N_phi):
    return np.array([(4*np.pi/N_phi)*i for i in range(N_phi)])
####

def get_psi(N_psi):
    return np.array([(4*np.pi/N_psi)*i for i in range(N_psi)])
####

def get_inner_prod(x, f1, f2, a, b):
    """
    inner product of 2 functions with weight function w(x)=1
    "x" is the symbol used in the symbolic expression of f1(x) and f2(x)
    """ 
    I = sp.Integral(f1*f2, (x, a, b))
    return I.doit()

def get_ws(N):
    """ Gaussian weights of order N """
    x = get_Legende_roots_x(N)
    w = [sp.symbols("w_{"+"{i}".format(i=i)+"}") for i in range(1, N+1)]
    LHS = sp.zeros(N+1, 1)
    for k in range(N):
        for s in range(N):
            LHS[k,0] += w[s] * sp.functions.special.polynomials.legendre(k, x[s])
    ####
    #
    RHS = sp.zeros(N+1, 1)
    xd = sp.symbols("x") # dummy variable
    p0 = sp.functions.special.polynomials.legendre(0, xd)
    inner_p0p0 = get_inner_prod(xd, p0, p0, -1, 1)
    RHS[0,0] = inner_p0p0
    #
    solution = sp.linsolve(LHS - RHS, w)

    # Extract the coefficients
    coefficients = solution.args[0]
    return list(coefficients)
##

def get_yi(alpha):
    """ returns yi of eq. (5) of https://arxiv.org/pdf/2304.02322.pdf from the angles """
    theta, phi, psi = alpha
    y0 = + np.cos(theta/2) * np.cos((phi+psi)/2)
    y1 = - np.cos(theta/2) * np.sin((phi+psi)/2)
    y2 = - np.sin(theta/2) * np.cos((-phi+psi)/2)
    y3 = + np.sin(theta/2) * np.sin((-phi+psi)/2)
    return [y0, y1, y2, y3]
####

def get_distance_S3(y1, y2):
    return np.linalg.norm(np.array(y1)-np.array(y2))
####

