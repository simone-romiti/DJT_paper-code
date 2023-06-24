# testing the locality of the operator whan the truncation approached infinity

import sys
sys.path.insert(0, '../')

import su2_DJT.operators.operators as operators
from   su2_DJT.DJT.DJT_matrix import *
import su2_DJT.S3_sphere.partition as partition

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle
import os

N_g = operators.N_g # number of generators in the algebra

q_vals = [1/2]
La_dict = {1: [], 2: [], 3: []}
abs_La_dict = {1: [], 2: [], 3: []}

# loading the data if present
data_path = "data.pickle"
if os.path.exists(data_path):
    with open(data_path, "rb") as file:
        data = pickle.load(file)
        q_vals = data["q_vals"]
        q_next = max(q_vals) + 1/2 # next value of q
        q_vals.append(q_next)

        La_dict = data["La_dict"]
        abs_La_dict = data["abs_La_dict"]
    ####
####


q_max = max(q_vals) # maximum number of q considered so far

print("Building the matrices for the L_a")
for q in [q_vals[-1]]:
    print("q =", q)
    DJT = get_DJT(q)
    DJT_dag = get_DJT_dag(DJT = DJT)
    #
    for a in range(1, N_g+1):
        La = operators.get_La(a=a, q=q, V=DJT, V_inv=DJT_dag)
        La_dict[a].append(La)
        abs_La = np.abs(La)/partition.get_asympt_d3alpha(q=q)
        abs_La_dict[a].append(abs_La)
    ####
####

file = "data.pickle"
with open("data.pickle", "wb") as file:
    data = {"q_vals" : q_vals,  "La_dict": La_dict, "abs_La_dict": abs_La_dict}
    pickle.dump(data, file)
    file.close()
####

print("Creating the animation")
out_dir = "./heatmaps/"
num_frames = len(q_vals)

for a in range(1, N_g+1):
    # list of matrices for the animation
    matrices = abs_La_dict[a]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Initialize the image object
    im = ax.imshow(matrices[0], cmap='hot', interpolation='none')

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
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    ####
    out_file = os.path.abspath(out_dir + "/L{a}-qmax{q_max}.html".format(a=a, q_max=q_max))
    animation.save(out_file, writer='html')    
####

