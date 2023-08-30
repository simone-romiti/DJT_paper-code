# 1-dimensional exact derivative for trigonometric polynomials of finite degree

# https://link.springer.com/content/pdf/10.1007/BF02748342.pdf?pdf=button

import numpy as np
import sympy as sp
from sympy.physics.quantum import TensorProduct as tp
from sympy.physics.quantum.spin import WignerD as WignerD
from sympy.physics.quantum.spin import Rotation as Wigner_rot

class D_alpha:
    def __init__(self, tau, x):
        self.tau = tau
        self.x = x
        self.n = len(list(self.x))
        y = sp.symbols("y")
        self.sprime = lambda x0: self.symb_s_prime(y).subs(y, x0).evalf()
    ####

    def symb_s(self, y: sp.Symbol):
        res = 1
        for xl in self.x:
            res *= sp.sin((y - xl)/2)
        ####
        return res
    ####
    def symb_s_prime(self, y: sp.Symbol):
        return sp.diff(self.symb_s(y), y)
    ####

    def get_nodes_matr(self):
        """ 
        diagonal matrix of the values of the nodes

        Returns: sympy Matrix
        """
        x = self.x
        n = self.n
        M = sp.zeros(n,n)
        for i in range(n):
            M[i,i] = x[i]
        ####
        return np.matrix(M)
    ####

    def get_Dtilde_angle(self):
        """
        Differential operator of eq. 3(a-b) of https://link.springer.com/content/pdf/10.1007/BF02748342.pdf?pdf=button

        Returns: symbolic sympy matrix
        """
        x = self.x
        n = len(list(x))
        D = sp.zeros(n, n)
        fact = sp.pi/self.tau
        for j in range(n):
            for k in range(n):
                if j==k:
                    s = 0
                    for m in range(n):
                        if m != j:
                            s += sp.cot(sp.pi*(x[j]-x[m])/self.tau)
                    ####
                    D[j, k] = s
                else:
                    D[j, k] = 1/sp.sin(sp.pi*(x[j]-x[k])/self.tau)
                ####
            ####
        ####
        D = fact*D
        return np.matrix(D)
    ####

    def get_matrix(self):
        x = self.x
        n = self.n
        D_alpha = np.zeros(shape = (n, n))
        for i in range(n):
            for j in range(n):
                res = 0
                if i==j:
                    for l in range(n):
                        if l != i:
                            res += 1 / np.tan((x[i] - x[l])/2)
                        ####
                    ####
                ####
                else:
                    res = (self.sprime(x[i])/self.sprime(x[j]))/np.sin((x[i]-x[j])/2)
                ####
                D_alpha[i,j] = res
            ####
        ####
        return (np.pi/self.tau) * np.matrix(D_alpha)
    ####
####


"""
Example:
N_pts = 15
x = np.array([i * 2*np.pi/N_pts for i in range(N_pts)])
f1 = np.matrix(np.sin(x)).transpose()
f2 = np.matrix(np.cos(x)).transpose()

D_obj = D_alpha(tau = 2*np.pi, x=x)
D = D_obj.get_matrix()

print("You should get the zero vector (up to machine precision)")
print(D*f1 - f2)
"""