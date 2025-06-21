import numpy as np


def H_ising(spin_list_matrices):
    """
    Retorna el hamiltoniano isin, el spin sz0 interactua via J con el resto de
    los espines
    spin_list_matrices: Lista de operadores de matrices de spin 
    Returns:  
    """
    h_ising = 0
    if len(spin_list_matrices) == 1:
        h_ising = np.zeros((2,2), dtype=complex)
    else:
        sz0 = spin_list_matrices[0]
        sz_list = spin_list_matrices[1:]
        h_ising = -4* sz0 @ np.sum(sz_list,axis=0)

    return h_ising

def H_disorder(o_sx, o_sz, disorder_x, disorder_z):
    """
    o_sx: List of spin matrices along x-direction
    o_sz: List of spin matrices along z-direction
    disorder_x: Disorder vector along x-direction
    disorder_z: Disorder vector along z-direction
    """
    n = len(disorder_x)
    suma = np.zeros_like(o_sz[0])
    for i in range(n):
        suma += disorder_x[i]*o_sx[i] + disorder_z[i]*o_sz[i]
    return suma


def H_kick(operators, pulse_angle, deviation_angle):
    """
    operators: List of spim matrices where the pulse is applied
    pulse_angle: Angle of the pulse applied
    deviation_angle: Deviation of the pulse angle is applied
    """
    return (pulse_angle - deviation_angle) * np.sum(operators, axis=0)


def H_drive(spin_list_matrices):
    """
    spin_list_matrices: Lista de matrices de espin
    """
    return np.sum(spin_list_matrices,axis=0)
