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
