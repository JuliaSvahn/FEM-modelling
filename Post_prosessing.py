import math
import numpy as np
import os
def postprocess(U, F, K, L, Trans, E):
    """
    Post-process the results of a finite element analysis.
    U: Displacement np.array
    F: Force np.array
    K: Stiffness np.array
    L: Length np.array
    Trans: Transformation matrix np.array
    E: Young's modulus np.array
    """
    # Calculate the reaction forces

    R = K @ U - F

    Strain = strain(U, L, Trans)


    # Calculate the stress in each element
    Stress = E * Strain
    return R, Stress, Strain

def strain(U, L, trans):
    """
    Calculate the strain in each element.
    U: Displacement np.array
    L: Length np.array
    trans: Transformation matrix np.array
    """
    strain = np.array([])
    for element in range(len(trans)):
        
        strain_element = (np.linalg.norm(U[trans[element][1]] - U[trans[element][0]])) / L[element]
    


        # print(U[trans[element][1]], U[trans[element][0]], L[element], strain_element)
        strain = np.append(strain, strain_element)

    """
    strain = np.array()
    """
    return strain