
import sys
sys.path.insert(0, '../')

import su2_DJT.operators.operators as operators
import numpy as np

q = 3

U = operators.get_U(q)
U_dag = operators.get_Udag(U)

Up = operators.color_prod_links(U_dag, U)
N_alpha = Up[0][0].shape[0]

print("""
  The following output should be the tensor product of the identity on the Hilbert space 
  with the identity on color space, 
  namely the matrix: Up = [\mathrm{1}, 0], [0, mathrm{1}]
  """)

decimals = 15
print("rounding decimals", decimals)

print("Up_11: ", end="")
print(np.array_equal(np.eye(N_alpha), Up[0][0].round(decimals=decimals)))

print("Up_12: ", end="")
print(np.array_equal(np.zeros(shape=(N_alpha, N_alpha)), Up[0][1].round(decimals=decimals)))

print("Up_21: ", end="")
print(np.array_equal(np.zeros(shape=(N_alpha, N_alpha)), Up[1][0].round(decimals=decimals)))

print("Up_22: ", end="")
print(np.array_equal(np.eye(N_alpha), Up[1][1].round(decimals=decimals)))

## print(U * U_dag)
