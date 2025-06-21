from scipy.linalg import expm, eig
import numpy as np


def U0(H_t, beta, T, t_step):
    """
    H_t: Funcion Hamiltoniano dependiente del tiempo
    T: Periodo
    t_step: Tamanho para discretizar el tiempo
    """
    steps = int(T/t_step)
    N = H_t(0,beta).shape[0]
    U0 = np.eye(N, dtype=complex)
    dt = t_step

    for k in range(steps):
        t = dt * k
        Ht = H_t(t,beta)
        U0 = expm(-1j * Ht * dt) @ U0

    return U0
