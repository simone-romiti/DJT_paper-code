# https://link.springer.com/content/pdf/10.1007/BF02748342.pdf?pdf=button

import sympy as sp
import numpy as np

from su2_DJT.fd_trig_poly import D_alpha


def L1(theta, phi, psi):
    D_theta = D_alpha(tau = np.pi, x=theta).get_matrix()
    a = - sp.cot(theta) * sp.cos(phi) * sp.diff(expr, phi)
    b = - sp.sin(phi) * sp.diff(expr, theta)
    c = + (sp.cos(phi) / sp.sin(theta)) * sp.diff(expr, psi)
    res = (-sp.I) * (a+b+c)
    return res
####

def L2(theta, phi, psi, expr):
    a = - sp.cot(theta) * sp.sin(phi) * sp.diff(expr, phi)
    b = + sp.cos(phi) * sp.diff(expr, theta)
    c = + (sp.sin(phi) / sp.sin(theta)) * sp.diff(expr, psi)
    res = (-sp.I) * (a+b+c)
    return res
####

def L3(theta, phi, psi, expr):
    return (-sp.I) * sp.diff(expr, phi)
####

def Lplus(theta, phi, psi, expr):
    return L1(theta, phi, psi, expr) + (sp.I) * L2(theta, phi, psi, expr)
####

def Lminus(theta, phi, psi, expr):
    Lp = Lplus(theta, phi, psi, expr)
    return sp.conjugate(Lp)
####


# ##########
# ##########
# ##########

# def get_Dtheta(theta):
#     return get_Dangle(n, sp.pi, theta)
# ####

# def get_phi_vals(n):
#     return get_nodes_vals(n, 0, 2*sp.pi)
# ####

# def get_phi_matr(n):
#     return get_nodes_matr(n, 0, 2*sp.pi)
# ####

# def get_Dphi(n):
#     return get_Dangle(n, 2*sp.pi, get_phi_vals(n))
# ####

# def get_psi_vals(n):
#     return get_nodes_vals(n, 0, 4*sp.pi)
# ####

# def get_psi_matr(n):
#     return get_nodes_matr(n, 0, 4*sp.pi)
# ####

# def get_Dpsi(n):
#     return get_Dangle(n, 4*sp.pi, get_psi_vals(n))

# def apply_diag(fun, M):
#     """
#     Apply function to diagonal matrix

#     Returns sympy matrix of elements fM_{ij} = fun(M_{ij})
#     """
#     n = M.shape[0]
#     fM = sp.zeros(n,n)
#     for i in range(n):
#         fM[i,i] = fun(M[i,i])
#     ####
#     return fM
# ####

# def one_over(x):
#     return 1/x
# ####

# def get_xiR_1(n):
#     theta = get_theta_matr(n)
#     psi = get_psi_matr(n)
#     Id_phi = sp.eye(n)
#     cot_theta = apply_diag(sp.cot, theta)
#     sin_theta = apply_diag(sp.sin, theta)
#     sin_theta_inv = apply_diag(one_over, sin_theta)
#     cos_psi = apply_diag(sp.cos, psi)
#     sin_psi = apply_diag(sp.sin, psi)
#     Dtheta = get_Dtheta(n)
#     Dphi = get_Dphi(n)
#     Dpsi = get_Dpsi(n)
#     c1 = tp(-cot_theta   , Id_phi, cos_psi*Dpsi)
#     c2 = tp(Dtheta       , Id_phi, -sin_psi    )
#     c3 = tp(sin_theta_inv, Dphi  , cos_psi     )
#     return c1+c2+c3
# ####

# def get_xiR_2(n):
#     theta = get_theta_matr(n)
#     psi = get_psi_matr(n)
#     Id_phi = sp.eye(n)
#     cot_theta = apply_diag(sp.cot, theta)
#     sin_theta = apply_diag(sp.sin, theta)
#     sin_theta_inv = apply_diag(one_over, sin_theta)
#     cos_psi = apply_diag(sp.cos, psi)
#     sin_psi = apply_diag(sp.sin, psi)
#     Dtheta = get_Dtheta(n)
#     Dphi = get_Dphi(n)
#     Dpsi = get_Dpsi(n)
#     c1 = tp(-cot_theta   , Id_phi, sin_psi*Dpsi)
#     c2 = tp(Dtheta       , Id_phi, cos_psi     )
#     c3 = tp(sin_theta_inv, Dphi  , cos_psi     )
#     return c1+c2+c3
# ####

# def get_xiR_3(n):
#     return tp(sp.eye(n), sp.eye(n), get_Dpsi(n))
# ####

# def get_R1(n):
#     return -sp.I*get_xiR_1(n)
# ####

# def get_R2(n):
#     return -sp.I*get_xiR_2(n)
# ####

# def get_R3(n):
#     return -sp.I*get_xiR_3(n)
# ####

# def get_Rplus(n):
#     R1, R2 = get_R1(n), get_R2(n)
#     return R1 + sp.I*R2
# ####

# def get_Rminus(n):
#     R1, R2 = get_R1(n), get_R2(n)
#     return R1 - sp.I*R2
# ####

# ###############################àà

# def get_xiL_1(n):
#     theta = get_theta_matr(n)
#     phi = get_phi_matr(n)
#     Id_psi = sp.eye(n)
#     cot_theta = apply_diag(sp.cot, theta)
#     sin_theta = apply_diag(sp.sin, theta)
#     sin_theta_inv = apply_diag(one_over, sin_theta)
#     cos_phi = apply_diag(sp.cos, phi)
#     sin_phi = apply_diag(sp.sin, phi)
#     Dtheta = get_Dtheta(n)
#     Dphi = get_Dphi(n)
#     Dpsi = get_Dpsi(n)
#     c1 = -tp(sin_theta_inv, cos_phi     , Dpsi)
#     c2 = tp(Dtheta        , sin_phi     , Id_psi)
#     c3 = tp(cot_theta     , cos_phi*Dphi, Id_psi)
#     return c1+c2+c3
# ####

# def get_xiL_2(n):
#     theta = get_theta_matr(n)
#     phi = get_phi_matr(n)
#     Id_psi = sp.eye(n)
#     cot_theta = apply_diag(sp.cot, theta)
#     sin_theta = apply_diag(sp.sin, theta)
#     sin_theta_inv = apply_diag(one_over, sin_theta)
#     cos_phi = apply_diag(sp.cos, phi)
#     sin_phi = apply_diag(sp.sin, phi)
#     Dtheta = get_Dtheta(n)
#     Dphi = get_Dphi(n)
#     Dpsi = get_Dpsi(n)
#     c1 = tp(sin_theta_inv, sin_phi       , Dpsi)
#     c2 = tp(Dtheta        , cos_phi      , Id_psi)
#     c3 = -tp(cot_theta     , sin_phi*Dphi, Id_psi)
#     return c1+c2+c3
# ####

# def get_xiL_3(n):
#     return tp(sp.eye(n), get_Dphi(n), sp.eye(n))
# ####

# def get_L1(n):
#     return sp.I*get_xiL_1(n)
# ####

# def get_L2(n):
#     return sp.I*get_xiL_2(n)
# ####

# def get_L3(n):
#     return sp.I*get_xiL_3(n)
# ####

# def get_Lplus(n):
#     L1, L2 = get_L1(n), get_L2(n)
#     return L1 + sp.I*L2
# ####

# def get_Lminus(n):
#     L1, L2 = get_L1(n), get_L2(n)
#     return L1 - sp.I*L2
# ####

# def get_Lsquared(n):
#     L1, L2, L3 = get_L1(n), get_L2(n), get_L3(n)
#     return L1*L1 + L2*L2 + L3*L3
# ####

# #####################
# #####################
# #####################





# def Wigner_d(j, m1, m2, beta):
#     """ 
#     Returns the wigner d functions as defined in 
#     https://docs.sympy.org/latest/modules/physics/quantum/spin.html#sympy.physics.quantum.spin.WignerD
#     """
#     return Wigner_rot.D(j=j, m=m1, mp=m2, alpha=0, beta=beta, gamma=0).doit()
# ####


# def WignerD_vec(n, j, m, mu):
#     """
#     Vector tensor product representation of eq. A.7 of 
#     https://link.springer.com/content/pdf/bbm:978-3-540-29082-7/1
#     """
#     theta, phi, psi = get_theta_vals(n), get_phi_vals(n), get_psi_vals(n)
#     d_fun = lambda t: Wigner_d(j=j, m1=m, m2=mu, beta=t) 
#     f1 = fun1_to_vec(d_fun, theta)
#     f2 = fun1_to_vec(sp.exp, -sp.I*m*phi)
#     f3 = fun1_to_vec(sp.exp, -sp.I*mu*psi)
#     return tp(f1, f2, f3)
# ####
