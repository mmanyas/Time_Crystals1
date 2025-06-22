import numpy as np
from static_functions.statics_functions import*
from hamiltonian_models.hamiltonian import*
from evolution_time.dynamics_functions import*
from observables.observaveis import*


def simulation_no_avg(omega,
                      theta,
                      theta0,
                      theta_eps,
                      beta,
                      delta_b,
                      delta,
                      T,
                      dt,
                      Nspin_vals,
                      T_steps
                      ):
    """
    omega: Frequency of the external drive term
    theta: Angulo del pulso aplicado
    theta0: Angulo de rotacion inicial del estado cuantico
    theta_eps: Desviacion del angulo theta del pulso aplicado
    beta: Amplitud of the extrenal drive field
    delta_b: Change of beta
    delta: Size to limit disorder intervals
    T: Periodo Normal de floquet, periodo del pulso aplicado
    dt: Ancho del tiempo discretizado en la construccion del operador unitario U
    Nspin_vals: Vector of spin number of system
    T_steps: Number of times applied of floquet periods
    """
    Nspin_vals=[2]
    for n in Nspin_vals:
        sx,_,sz = create_spin_xyz_operators(n)
        hx      = disorder_x(n, delta)
        hz      = disorder_z(n, delta)
        H       = hamiltoniano(n, omega, hx, hz, sx, sz)
        Hkick   = H_kick(sx,theta,theta_eps)
        psi0    = spin_rotator(sx, theta0) @ initial_state_up(n)
        #UF      = U_floquet1T(H,Hkick,beta,T,dt)
        UF      = U_floquet2T(H,Hkick,beta,T,dt)
        UF2     = U_floquet2T(H,Hkick,beta+delta_b,T,dt)
        print("hx", hx)
        print()
        print("hz", hx)

        filename = f"resultN={n}.dat"
        fout = open(filename, "w")
        fout.write("# t\tmagnetization\tinfo_fisher\tmetodo_momentos\n")

        for l in range(1,T_steps+1):
            psit  = psi_expandido(UF,psi0,l)
            psit2 = psi_expandido(UF2,psi0,l)

            mz      = magnetization(psit,sz[0])
            fish    = 0
            FI_mom  = 0
            fout.write(f"{l}\t{mz:.16f}\t{fish.real:.16f}\t{FI_mom:.16f}\n")

        fout.close()


