import numpy as np


def magnetization(psi,Sz):
    """
    psi: Estado cuantico del sistema
    Sz: Componente z de la matriz de espines
    """
    return np.vdot(psi,Sz @ psi).real

