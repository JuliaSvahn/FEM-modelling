import math
import numpy as np
import os
def local_ESM(E, A, L):
    """
    Calculate the local coordinate system element stiffness matrix for a 2D truss element.

    Parameters:
    E : float
        Young's modulus
    A : float
        Cross-sectional area
    L : float
        Length

    Returns:
    K_local : numpy.ndarray [2x2]
        local stiffness matrix of element
    k_local_4x4 : numpy.ndarray [4x4]
        local stiffness matrix of element in 4x4 format for global assembly
    """
    
    #Guard to make sure we don't accidentally divide by 0
    if L <= 0:
        raise ValueError(f"Element length must be positive, got L = {L}")

    # Calculate the stiffness coefficient
    k = (E * A) / L



    # Define the local stiffness matrix for a 2D truss element
    K_local = np.array([[ k, -k],
                        [-k,  k]])

    k_local_4x4 = np.zeros((4, 4))
    k_local_4x4[0, 0] = K_local[0, 0]
    k_local_4x4[0, 2] = K_local[0, 1]
    k_local_4x4[2, 0] = K_local[1, 0]
    k_local_4x4[2, 2] = K_local[1, 1]

    return K_local, k_local_4x4
