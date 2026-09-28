import math
import numpy as np
import os
def Rotation(x1, x2, y1, y2):

    """
    Calculate the rotation matrix for a 2D truss element.
    inputs:
    x1, y1: coordinates of the first node
    x2, y2: coordinates of the second node
    
    returns:
    rotation_matrix: 4x4 numpy array representing the rotation matrix
    theta: angle of rotation in radians
    """

    dy = y2-y1
    dx = x2-x1
    theta = np.arctan2(dy, dx)

    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta), 0, 0], 
                                      [np.sin(theta), np.cos(theta), 0, 0],
                                      [0, 0, np.cos(theta), -np.sin(theta)],
                                      [0, 0, np.sin(theta), np.cos(theta)]])

    return rotation_matrix, theta