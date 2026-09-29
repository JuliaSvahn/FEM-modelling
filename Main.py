import math
import numpy as np
import os

import ESM
import GSM
import Rotation_matrix
import Solve_system
import Post_prosessing



nodes = {}
elements = {}
BCs = {}
loads = {}
info = ""

filepath = "model 1.txt"
with open(filepath, "r") as f:

    for line in f:
        newline = line.strip()
        if not newline:
            continue
        if newline[0] == "[":
            info = newline[1:-1].strip()
            continue
        
        column = []
        for item in newline.split(","):
            column.append(item.strip())
        
        if info == "NODES":
            xcoord = float(column[1])
            ycoord = float(column[2])
            
            nodes[int(column[0])] = (xcoord, ycoord)
        
        elif info == "ELEMENTS":
            node1 = int(column[1])
            node2 = int(column[2])
            E = float(column[3])
            poisson = float(column[4])
            area = float(column[5])
            
            elements[int(column[0])] = (node1, node2, E, poisson, area)
        
        elif info == "BCS":
            BCx = int(column[1])
            BCy = int(column[2])
            
            BCs[int(column[0])] = (BCx, BCy)
        
        elif info == "LOADS":
            forcex = float(column[1])
            forcey = float(column[2])
            
            loads[int(column[0])] = (forcex, forcey)

# calculations

# ????? should the indexes be adjusted here???

# getting the length for each element then using element stiffness matrix
element_lengths = {}
element_angles = {}
rotation_matrices = {}

n_dof = 2 * len(nodes)

# now that we are done with bullshit lists its time for the real programming
# Convert each dict's values into a 2D NumPy array
nodes_arr = np.array(list(nodes.values()))
elements_table = np.array(list(elements.values()))
elements_arr = elements_table[:, :2].astype(int) - 1
E_arr = elements_table[:, 2]
A_arr = elements_table[:, 4]

k_global = GSM.assemble_truss_stiffness(nodes_arr, elements_arr, E_arr, A_arr, dof_node = 2)

# load vector
F = np.zeros(n_dof)
for nodeid, (forcex, forcey) in loads.items():
    i = nodeid - 1
    F[2 * i] += forcex # i don't really understand this section
    F[2 * i + 1] += forcey

fixedBC = []
for nodeid, (BCx, BCy) in BCs.items():
    i = nodeid - 1
    if BCx == 1:
        fixedBC.append(2 * i)
    if BCy == 1:
        fixedBC.append( 2 * i + 1)

u, RF = Solve_system.solve_gsm(k_global, F, fixedBC)

L_arr = []
T_list = []
for i, j in elements_arr:
    x1, y1 = nodes_arr[i]
    x2, y2 = nodes_arr[j]
    L_arr.append(np.hypot(x2 - x1, y2 - y1))
    T, theta = Rotation_matrix.Rotation(x1, x2, y1, y2)
    T_list.append(T)

R, stress, strain = Post_prosessing.postprocess(u, F, k_global, L_arr, T_list, E_arr, elements_arr)

print("Displacement: ", u)
print("Reaction forces: ", R)
print("Stresses: ", stress)
print("Strains: ", strain)
