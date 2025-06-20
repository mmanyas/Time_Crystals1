import numpy as np

sigma_x = np.array([[0,1],[1,0]], dtype=complex)/2
sigma_y = np.array([[0,-1j],[1j,0]], dtype=complex)/2
sigma_z = np.array([[1,0],[0,-1]], dtype=complex)/2

identity = np.eye(2, dtype=complex)

def sigma_i(direction):

    result = identity
    match direction:
        case 'x':
            result = sigma_x
        case 'y':
            result = sigma_y
        case 'z':
            result = sigma_z

    return result



def build_pmatrix(N, direction):
    """
    Devuelve la matriz de Pauli direccionada a un sitio específico en una cadena de N espines.

    Args:
        N (int): Número de espines
        site (int): Sitio donde se aplica el operador (0-indexado)
        direction (str): 'x', 'y' o 'z'

    Returns:
        np.ndarray: Matriz 2^N x 2^N del operador deseado
    """
    if direction == 'x':
        op = sigma_x
    elif direction == 'y':
        op = sigma_y
    elif direction == 'z':
        op = sigma_z
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








def create_Sxyz(n):
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
                sx = np.kron(sx, sigma_x) # tensor product
                sy = np.kron(sy, sigma_y)
                sz = np.kron(sz, sigma_z)
            else:
                sx = np.kron(sx, identity)
                sy = np.kron(sy, identity)
                sz = np.kron(sz, identity)

        sx_list.append(sx)
        sy_list.append(sy)
        sz_list.append(sz)

    return sx_list, sy_list, sz_list

