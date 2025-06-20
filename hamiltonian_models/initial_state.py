import numpy as np
from scipy.linalg import expm

def build_initial_state(N, spin_up_down,operators,initial_rotation):
    """
    Creates an initial quantum state by rotating a product state of N qubits.
    
    Parameters:
    N (int): Number of qubits
    operators (list of ndarray): List of operator matrices for each qubit
    initial_rotation (float): Rotation angle on initial state
    
    Returns:  ndarray: The initial quantum state vector
    """

    total_operator = -1j * initial_rotation * sum(operators)
    rotate         = expm(total_operator)
    
    # Create initial state |0⟩⊗|0⟩⊗...⊗|0⟩ (N times)
    up   = np.array([[1],[0]])
    down = np.array([[0],[1]])
    psi0 = []

    for _ in range(N):
        if spin_up_down=='up':
            psi0 = np.kron(up, np.array([[1],[0]]))  # |0⟩ = [[1], [0]]
        elif spin_up_down=='down':
            psi0 = np.kron(down, np.array([[0],[1]]))
        else:
            print("Choose the spin state: 'up' or 'down'")
    
    return rotate @ psi0
