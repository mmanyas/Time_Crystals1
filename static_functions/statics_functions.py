import numpy as np

sigma_x = np.array([[0,1],[1,0]], dtype=complex)
sigma_y = np.array([[0,-1j],[1j,0]], dtype=complex)
sigma_z = np.array([[1,0],[0,-1]], dtype=complex)

identity = np.eye(2, dtype=complex)

def sigma_i(direction):
    """
    direction: axis direction of the Pauli matrix, can be x,y,z
    return: return Pauli spin matrix
    """

    result = identity
    match direction:
        case 'x':
            result = sigma_x
        case 'y':
            result = sigma_y
        case 'z':
            result = sigma_z

    return result



def create_spin_operators(N, direction):
    """
    Devuelve la matriz de Pauli direccionada a un sitio específico en una cadena de N espines.
    Args:
        N (int): Número de espines del sistema
        direction (str): 'x', 'y' o 'z'
        Returns: list of np.ndarray, Matriz 2^N x 2^N del operador deseado
    """
    
    if direction == 'x':
        op = sigma_x/2
    elif direction == 'y':
        op = sigma_y/2
    elif direction == 'z':
        op = sigma_z/2
    else:
        raise ValueError("Dirección inválida. Usa 'x', 'y' o 'z'.")

    result = []

    for i in range(N):
        matriz = 1

        for j in range(N):
            if j==i:
                matriz = np.kron(matriz, op)
            else:
                matriz = np.kron(matriz, identity)

        result.append(matriz)

    
    return result





def create_spin_xyz_operators(n):
    """Crea matrices de spin para un sistema de n partículas.
    Args:
        n: Número de partículas/spins en el sistema.

    Returns:
        Tres listas conteniendo las matrices Sx, Sy, Sz para cada partícula,
        donde cada matriz está en el espacio de Hilbert completo del sistema.
    """
    sx_list = []
    sy_list = []
    sz_list = []

    for i in range(n):
        sx = 1
        sy = 1
        sz = 1

        for j in range(n):
            if j == i:
                sx = np.kron(sx, sigma_x/2) # tensor product
                sy = np.kron(sy, sigma_y/2)
                sz = np.kron(sz, sigma_z/2)
            else:
                sx = np.kron(sx, identity)
                sy = np.kron(sy, identity)
                sz = np.kron(sz, identity)

        sx_list.append(sx)
        sy_list.append(sy)
        sz_list.append(sz)

    return sx_list, sy_list, sz_list

