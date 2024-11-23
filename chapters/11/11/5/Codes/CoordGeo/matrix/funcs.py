#Code by GVV Sharma
#July 26, 2024
#released under GNU GPL
#Adapted from Chatgpt
#Matrix functions

import numpy as np
import numpy.linalg as LA

def polyvalm(f, A):
    """
    Evaluate a matrix polynomial at matrix A.
    
    Parameters:
    f : numpy.ndarray
        Coefficients of the polynomial in decreasing order of power.
    A : numpy.ndarray
        The matrix at which the polynomial is evaluated.
    
    Returns:
    numpy.ndarray
        The result of evaluating the polynomial at matrix A.
    """
    
    # Initialize the result matrix as zero matrix
    result = np.zeros_like(A, dtype=np.float64)
    
    # Evaluate the polynomial by summing up each term
    for i, coeff in enumerate(f):
        result += coeff * LA.matrix_power(A, len(f) - 1 - i)
    
    return result

