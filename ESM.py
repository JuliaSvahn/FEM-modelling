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
    """
    
    #Guard to make sure we don't accidentally divide by 0
    if L <= 0:
        raise ValueError(f"Element length must be positive, got L = {L}")

    # Calculate the stiffness coefficient
    k = (E * A) / L

    # Define the local stiffness matrix for a 2D truss element
    K_local = np.array([[ k, -k],
                        [-k,  k]])

    return K_local
