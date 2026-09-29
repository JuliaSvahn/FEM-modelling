import numpy as np
def postprocess(U, F, K, L, Trans, E, elements):
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

    Strain = strain(U, L, Trans, elements)


    # Calculate the stress in each element
    Stress = E * Strain
    return R, Stress, Strain


def strain(U, L, trans, elements):
    """
    Calculate the strain in each element.
    U: Displacement np.array
    L: Length np.array
    trans: Transformation matrix np.array
    elements: Element connectivity array
    """
    eps = np.zeros(len(elements))
    for e, (i, j) in enumerate(elements):
        dofs = [2*i, 2*i+1, 2*j, 2*j+1]   
        u_local = trans[e].T @ U[dofs]           
        eps[e] = (u_local[2] - u_local[0]) / L[e]  
    return eps
