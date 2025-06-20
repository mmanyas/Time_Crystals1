import numpy as np
from scipy.linalg import expm 

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
    Args:
        N (int): Número de espines del sistema
        direction (str): 'x', 'y' o 'z'
        Returns: list of np.ndarray, matrix 2^N x 2^N del operador deseado
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
    """ Crea matrices de spin para un sistema de n partículas.
    Args:
        n: Número de partículas/spins en el sistema.
    Returns: Sx_list, Sy_list, Sz_list. if n=2, sx_list = [sx0, sx1] ....
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

def disorder_x(Nl, deltal):
    """
    create a disorder vector in x-direction between [-delta/2:delta/2]
    Nl: Number of sipins of the system
    deltal: Limits the values of the randomly generated disorder vector
    Returns: Returns a disorder vector of Nl size
    """
    return np.random.uniform(-deltal/2, deltal/2, Nl)


def disorder_y(Nl, deltal):
    """
    create a disorder vector in y-direction between [-delta/2:delta/2]
    Nl: Number of sipins of the system
    deltal: Limits the values of the randomly generated disorder vector
    Returns: Returns a disorder vector of Nl size
    """
    return np.random.uniform(-deltal/2, deltal/2, Nl)


def disorder_z(Nl, deltal):
    """
    create a disorder vector in z-direction between [-delta/2:delta/2]
    Nl: Number of sipins of the system
    deltal: Limits the values of the randomly generated disorder vector
    Returns: Returns a disorder vector of Nl size
    """
    return np.random.uniform(0, deltal/2, Nl)


def initial_state_up(Nl):
    """
    Create initial quantum state with all spin up of Nl spin-system
    Nl: Number of spins 
    Returns: Returns a column vector 
    """
    up   = np.array([[1],[0]], dtype='complex')
    psi0 = 1
    for _ in range(Nl):
        psi0 = np.kron(psi0,up)
    return psi0


def initial_state_down(Nl):
    """
    Create initial quantum state with all spin down of Nl spin-system
    Nl: Number of spins 
    Returns: Returns a column vector 
    """
    down = np.array([[0],[1]], dtype='complex')
    psi0 = 1
    for _ in range(Nl):
        psi0 = np.kron(psi0,down)
    return psi0


def spin_rotator(opr, angle):
    """
    This function create a spin operador rotator
    opr: Matrix rotation, can be Pauli matrices
    angle: Angle of rotations
    """
    return expm( -1j * angle * sum(opr) )


def check_rotation(rot):
    inver = rot.conj().T
    check = inver @ rot
    N = rot.shape[0]

    is_unitary = np.allclose(check, np.eye(N))
    print("Is the rotation unitary? ", is_unitary)


