import math
import numpy as np
import os
def solve_gsm(K, loads, fixed_dofs):
    """
    Solving GSM system for displacements and reactions.

    Elements:
    K -> Global stiffness matrix
    loads -> Load vector applied
    fixed_dofs -> dofs used to solve system
    """
    deg_node = 2 # degrees of freedom per node
    n_dof = K.shape[0] #No of dofs in total

    F = np.zeros(n_dof)
    for load in loads:
        dof_index = int(load[0]) - 1  # Convert to zero-based index
        force_value_x = load[1]
        force_value_y = load[2] # obtain the force values from the load array
        F[dof_index * deg_node] = force_value_x # assign the components to right indices as 2 dof means times 2
        F[dof_index * deg_node + 1] = force_value_y
        
    fixed_dofs = np.asarray(fixed_dofs) #No of dofs for solving the system
    free_dofs = np.setdiff1d(np.arange(n_dof), fixed_dofs) #All dofs that aren't required for the 

    u = np.zeros(n_dof)
    u[free_dofs] = np.linalg.solve(K[np.ix_(free_dofs, free_dofs)], F[free_dofs])

    RF = K @ u - F     
    return u, RF