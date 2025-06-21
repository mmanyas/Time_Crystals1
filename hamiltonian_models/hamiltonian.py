import numpy as np


def H_ising(spin_list_matrices):
    """
    Retorna el hamiltoniano isin, el spin sz0 interactua via J con el resto de
    los espines
    spin_list_matrices: Lista de operadores de matrices de spin 
    Returns:  
    """
    sz0 = spin_list_matrices[0]
    sz_list = spin_list_matrices[1:]

    return -4* sz0 @ np.sum(sz_list,axis=0)


def H_disorder(o_sx, o_sz, disorder_x, disorder_z):
    """
    o_sx: List of spin matrices along x-direction
    o_sz: List of spin matrices along z-direction
    disorder_x: Disorder vector along x-direction
    disorder_z: Disorder vector along z-direction\
    """
    n = len(disorder_x)
    suma = np.zeros_like(o_sz[0])
    for i in range(n):
        suma += disorder_x[i]*o_sx[i] + disorder_z[i]*o_sz[i]
    return suma
