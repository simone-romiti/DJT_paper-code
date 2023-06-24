# testing the locality of the operator whan the truncation approached infinity

import operators
from DJT_matrix import *
import partition

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle


q_vals = [k/2 for k in range(1, 8)]
q_max = max(q_vals)
N_g = operators.N_g # number of generators in the algebra

abs_La_dict = {1: [], 2: [], 3: []}

print("Building the matrices for the L_a")
for q in q_vals:
    print("q =", q)
    DJT = get_DJT(q)
    DJT_dag = get_DJT_dag(DJT = DJT)
    #
    for a in range(1, N_g+1):
        La = operators.get_La(a=a, q=q, V=DJT, V_inv=DJT_dag)
        abs_La = np.abs(La)/partition.get_asympt_d3alpha(q=q)
        abs_La_dict[a].append(abs_La)
    ####
####

with open('abs_La.pickle', 'wb') as file:
    pickle.dump(abs_La_dict, file)

print("Creating the animation")
num_frames = len(q_vals)
for a in range(1, N_g+1):
    # list of matrices for the animation
    matrices = abs_La_dict[a]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Initialize the image object
    im = ax.imshow(matrices[0], cmap='hot', interpolation='nearest')

    cbar = fig.colorbar(im)
    # Define the update function for animation
    def update(frame):
        im.set_array(matrices[frame])
        ax.set_title("q = {q}".format(q = q_vals[frame]))
        return im,
    ####

    # Create the animation
    animation = FuncAnimation(fig, update, frames=num_frames, interval=500, blit=True)

    # Save the animation as HTML
    animation.save("L{a}-qmax{q_max}.html".format(a=a, q_max=q_max), writer='html')    
####

