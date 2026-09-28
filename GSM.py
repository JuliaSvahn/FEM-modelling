import math
import numpy as np
import os

def transform_to_global(K_local, c, s):
    """
    Rotate a 2x2 local (axial) element stiffness matrix into the
    4x4 global stiffness matrix of a 2D truss element.
    inputs:
        K_local : 2x2 local stiffness matrix
        c, s    : cos(theta), sin(theta) of the element orientation
    """
    T = np.array([[c, s, 0, 0],
                  [0, 0, c, s]])
    return T.T @ K_local @ T



def assemble_truss_stiffness(nodes, elements, E, A, dof_node=2):
    """
    nodes    : array-like (n_nodes, 2)        -> [x, y] per node
    elements : array-like (n_elements, 2) int -> 0-indexed [node_i, node_j]
    E, A     : array-like (n_elements,)       -> per-element properties
    """
    K_global_list = []

    for e, (i, j) in enumerate(elements):
        xi, yi = nodes[i]
        xj, yj = nodes[j]
        dx, dy = xj - xi, yj - yi
        L = np.hypot(dx, dy)
        c, s = dx / L, dy / L

        K_local = local_ESM(E[e], A[e], L)                # 2x2, local axial
        K_global_e = transform_to_global(K_local, c, s)   # 4x4, global x/y
        K_global_list.append(K_global_e)

    K = GSM_assembly(elements, nodes, dof_node, K_global_list)
    return np.asarray(K, dtype=float)

def GSM_assembly(elements, nodes, K_local, dof_node):
    """
    Function assembling Global Stiffness Matrix (GSM) for a 2D truss structure.
    From other codes we have the element stiffness matrices, this code assembles them.
    """
    n_dof = dof_node * len(nodes)  # Total degrees of freedom + no of col's and rows in GSM

    K = np.zeros((n_dof, n_dof))  # Initialize the global stiffness matrix (square)

    for e, (i, j) in enumerate(elements):
        # Global DOF numbers of this element: node n owns DOFs n*dof_node ... n*dof_node + dof_node - 1
        dofs = np.concatenate([np.arange(i * dof_node, (i + 1) * dof_node),
                               np.arange(j * dof_node, (j + 1) * dof_node)])

        # Scatter-add the element matrix into the right rows/columns
        K[np.ix_(dofs, dofs)] += K_local[e]

    return K
