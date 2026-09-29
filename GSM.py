import math
import numpy as np
import os

def assemble_truss_stiffness(nodes, elements, E, A, dof_node=2):
    """
    nodes    : array-like (n_nodes, 2)        -> [x, y] per node
    elements : array-like (n_elements, 2) int -> 0-indexed [node_i, node_j]
    E, A     : array-like (n_elements,)       -> per-element properties
    """
    n_dof = dof_node * len(nodes)
    K = np.zeros((n_dof, n_dof)) # global stiffness matrix

    for e, (i, j) in enumerate(elements):
        xi, yi = nodes[i]
        xj, yj = nodes[j]
        dx, dy = xj - xi, yj - yi
        L = np.hypot(dx, dy)

        K_local, k4 = local_ESM(E[e], A[e], L)
        R, theta = Rotation (xi, xj, yi, yj)
        K_global = R @ k4 @ R.T  #switch to global coordinates

        dofs = np.concatenate([np.arange(i * dof_node, (i + 1) * dof_node),
                               np.arange(j * dof_node, (j + 1) * dof_node)])

        for a in range(len(dofs)):
            for b in range(len(dofs)):
                K[dofs[a], dofs[b]] += K_global[a, b] # add each element's stiffness matrix to the global matrix

    return K
