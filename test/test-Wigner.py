# Test of the orthogonality of Wigner D and d functions

import sympy as sp
from   sympy.physics.quantum.spin import WignerD

import sys
sys.path.insert(0, '../')

from su2_DJT.S3_sphere.partition import *

def Wignerd(j, m, mu, theta):
    return WignerD(j, m, mu, 0.0, theta, 0.0).doit()
####

theta = sp.symbols("theta")
q = 20
check_truncation(q)
N = int(2*q + 1)
theta_S3 = get_theta(N)
w = get_ws(N)
m, mu = 5/2, -17/2
for j1 in [l1/2 for l1 in range(0, q)]:
    if abs(m) > j1:
        continue
    ####
    for j2 in  [l2/2 for l2 in range(0, q)]:
        if abs(mu) > j2:
            continue
        ####
        b1, b2 = is_half_integer(j1), is_half_integer(j2)
        b3, b4 = is_half_integer(m), is_half_integer(mu)
        B1 = (b1 and b2 and b3 and b4)
        B2 = not (b1 or b2 or b3 or b4)
        if not (B1 or B2):
            continue # j1 and j2 have to be both integers or half integers
        ####
        f = sp.sin(theta) 
        f = f * Wignerd(sp.Rational(j1), sp.Rational(m), sp.Rational(mu), theta) 
        f = f * Wignerd(sp.Rational(j2), sp.Rational(m), sp.Rational(mu), theta) 
        I = sp.Integral(f, (theta, 0, sp.pi))
        I = (j1 + 1/2) * (I.doit().evalf(14, chop=True))
        print("Integral:", sp.Rational(j1), sp.Rational(j2), I)
        tot = 0.0
        for s in range(N):
            part_s = w[s]
            part_s = part_s * Wignerd(sp.Rational(j1), sp.Rational(m), sp.Rational(mu), theta_S3[s]) 
            part_s = part_s * Wignerd(sp.Rational(j2), sp.Rational(m), sp.Rational(mu), theta_S3[s])
            tot += part_s
        ####
        tot = (tot*(j1 + 1/2)).evalf(14, chop = True)
        print("Sum:", sp.Rational(j1), sp.Rational(j2), tot)
    ####
####


