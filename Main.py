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

filepath = "model 1.txt"
with open(filepath, "r") as f:
    
    info = ""
    
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
            nodeid = int(column[0])
            xcoord = float(column[1])
            ycoord = float(column[2])
            
            nodes[nodeid] = nodeid, xcoord, ycoord
        
        elif info == "ELEMENTS":
            elid = int(column[0])
            node1 = int(column[1])
            node2 = int(column[2])
            E = float(column[3])
            poisson = float(column[4])
            area = float(column[5])
            
            elements[elid] = elid, node1, node2, E, poisson, area
        
        elif info == "BCS":
            BCid = int(column[0])
            BCx = int(column[1])
            BCy = int(column[2])
            
            BCs[BCid] = BCid, BCx, BCy
        
        elif info == "LOADS":
            loadid = int(column[0])
            forcex = float(column[1])
            forcey = float(column[2])
            
            loads[loadid] = loadid, forcex, forcey

# calculations

# getting the length for each element then using element stiffness matrix
element_lengths = {}
element_angles = {}
rotation_matrices = {}
k_local_list = []


for elid, element in elements.items():

    elid, node1, node2, E, poisson, area = element

    x1, y1 = nodes[node1][1], nodes[node1][2]
    x2, y2 = nodes[node2][1], nodes[node2][2]

    dx = x2 - x1
    dy = y2 - y1

    L = np.sqrt(dx**2 + dy**2)

    element_lengths[elid] = L

    K_local, k_local_4x4 = ESM.local_ESM(E, area, L)
    
    T, theta = Rotation_matrix.Rotation(x1, x2, y1, y2)
    rotation_matrices[elid] = T
    element_angles[elid] = theta

    k_local_list.append(T @ k_local_4x4 @ T.T)

# now that we are done with bullshit lists its time for the real programming
# Convert each dict's values into a 2D NumPy array
nodes = np.array(list(nodes.values()))
elements = np.array(list(elements.values()))
BCs = np.array(list(BCs.values()))
loads = np.array(list(loads.values()))

k_global = GSM.GSM_assembly(elements, nodes, k_local_list, dof_node = 2) #hello
