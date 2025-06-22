from scipy.linalg import expm, eig
import numpy as np


def U_floquet1T(H_t, H_kick, beta, T, t_step):
    """
    H_t: Funcion Hamiltoniano dependiente del tiempo
    T: Periodo
    t_step: Tamanho para discretizar el tiempo
    """
    steps = int(T/t_step)
    N = H_t(0,beta).shape[0]
    U0 = np.eye(N, dtype=complex)
    dt = t_step

    for k in range(steps+1):
        t = dt * k
        Ht = H_t(t,beta)
        if t==T:
            Ht += H_kick
            print("entro cuando t = T = ", T)
        U0 = expm(-1j * Ht * dt) @ U0

    return U0


def U_floquet2T(H_t, H_kick, beta, T, t_step):
    """
    H_t: Funcion Hamiltoniano dependiente del tiempo
    T: Periodo
    t_step: Tamanho para discretizar el tiempo
    """
    steps = int(T/t_step)
    N = H_t(0,beta).shape[0]
    U0 = np.eye(N, dtype=complex)
    dt = t_step

    for k in range(2*steps+1):
        t = dt * k
        Ht = H_t(t,beta)
        if t==T and t==2*T:
            Ht += H_kick
            print("entro cuando t = T = ", T)
        U0 = expm(-1j * Ht * dt) @ U0

    return U0


def psi_expandido(U, psi0, n):
    """
    U: Operador de evolucion unitario a ser expandido
    psi0: Estado inicial del sistema
    n: Numero de pulsos de floquet aplicados
    """
    evals, evecs = eig(U)
    coef = [ np.vdot(evecs[:,i],psi0) for i in range(evecs.shape[1])]
    psit = np.zeros_like(psi0)
    for i in range(len(coef)):
        psit += evals[i]**n * coef[i] * evecs[:,i].reshape(-1,1)

    return psit
