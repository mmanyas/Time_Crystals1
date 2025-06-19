import numpy as np
from scipy.linalg import expm,eig


def create_spin_operators(N, direction):
    """
    Crea una lista de operadores de espín para N partículas en una dirección específica.
    
    Args:
        N (int): Número de espines/partículas
        direction (str): Dirección del espín ('x', 'y' o 'z')
    
    Returns:
        list: Lista de operadores de espín (matrices numpy) para cada partícula en la dirección especificada
    """
    # Matrices de Pauli
    pauli_matrices = {
        'x': np.array([[0, 1], [1, 0]]) / 2,
        'y': np.array([[0, -1j], [1j, 0]]) / 2,
        'z': np.array([[1, 0], [0, -1]]) / 2
    }
    
    I = np.eye(2)				# Matriz identidad
    
    if direction.lower() not in pauli_matrices:	# validar la direccion 
        raise ValueError("Dirección debe ser 'x', 'y' o 'z'")
    
    sigma = pauli_matrices[direction.lower()]
    spin_list = []
    
    for i in range(N):
        operator = [sigma if j == i else I for j in range(N)]
        current_op = 1
  
        for op in operator:
            current_op = np.kron(current_op, op)
        spin_list.append(current_op)
    
    return spin_list
    

#===============================================================================================    
#==========================================Hamiltoniano=========================================
def Hamiltonian(time_function, J, N, random_list, terms):
    """ 
    Retorna una función H(t) evaluable y una lista H_list compatible con QuTiP. 
    time_function: función f(t, beta, omega)
    J: acoplamiento de Ising
    N: número de espines
    delta: campo local aleatorio
    terms: tipo de sistema (0: Ising, 1: +campo, 2: +drive, 3: Ising+drive)
    """

    sx = create_spin_operators(N, 'x')
    sy = create_spin_operators(N, 'y')
    sz = create_spin_operators(N, 'z')    

#    np.random.seed(60)
 #   hx = np.random.uniform(-delta/2, delta/2, N)
 #   hz = np.random.uniform(0, delta, N)
    hx = random_list[0]
    hz = random_list[1]

    H_ising = -J * sz[0] * sum(sz[1:])
    H_pert  = sum(hx[i] * sx[i] + hz[i] * sz[i] for i in range(N))
    H_drive = sum(sz)

    # H_list para qutip
    if terms == 0:
        H_list = [H_ising]
    elif terms == 1:
        H_list = [H_ising + H_pert]
    elif terms == 2:
        H_list = [H_ising + H_pert, [H_drive, lambda t, args: time_function(t, args['beta'], args['omega'])]]
    elif terms == 3:
        H_list = [H_ising, [H_drive, lambda t, args: time_function(t, args['beta'], args['omega'])]]
    else:
        raise ValueError("El parámetro 'terms' debe ser 0, 1, 2 o 3.")

        
        
    # Función evaluable H(t)
    def H_func(t, beta, omega):
        if terms == 0:
            print("termssssssssssssssssssssssssssss = ", terms, "beta = ", beta, "omega",omega)
            return H_ising
        elif terms == 1:
            return H_ising + H_pert
            print("termssssssssssssssssssssssssssss = ", terms, "beta = ", beta, "omega",omega)
            print("function = ", time_function(t,beta,omega))
        elif terms == 2:
            print("termssssssssssssssssssssssssssss = ", terms, "beta = ", beta, "omega",omega)
            print("function = ", time_function(t,beta,omega))
            return H_ising + H_pert + time_function(t, beta, omega) * H_drive
        elif terms == 3:
            return H_ising + time_function(t, beta, omega) * H_drive
            print("termssssssssssssssssssssssssssss = ", terms, "beta = ", beta, "omega",omega)
        
        
    return H_func

#=========================================================================================
#========================== operador de evolucion unitario ===============================

def build_U2(hamiltonianl, Jl, Nl, termsl, Tl, dtl, betal, omegal, thetal, operatorsl):
		#hamiltonianl,Jl,Nl,H_termsl, T,dtl,betal,omegal,sx
    
    #H_func, H_list = Hamiltonian2(f_time,J,Ns,delta,terms)
    stepsl = int(Tl/dtl)
    times = np.linspace(0,Tl,stepsl)

    U_t = np.eye(2**Nl,dtype=complex)
    #print("betal = ",  betal)
    #print()
    for tk in range(2*stepsl + 1):
        t = tk*dt
        H = hamiltonianl(t,betal,omegal)
        U_t = expm( -1j * H * dtl ) @ U_t
        #print(H)        
        if t==T or t==2*T:
            U_t = expm(-1j * thetal * sum(operatorsl)) @ U_t
            
    #U_t = U_t, dims=[[2]*Ns,[2]*Ns] )
    
    return U_t
    
#========================================================================
#========================== Sz spin basis ===============================    
def szbasis(Ns, direction):
    """
    Genera el estado producto de N espines (qubits) todos en |↑⟩ o |↓⟩.
    
    Args:
        direction (str): "up" para |↑⟩ o "down" para |↓⟩.
        N (int): Número de espines.
        
    Returns:
        np.ndarray: Vector columna de shape (2^N, 1) con el estado.
    """
    
    if direction == "up":
        single_spin = np.array([[1], [0]])  # |↑⟩
    elif direction == "down":
        single_spin = np.array([[0], [1]])  # |↓⟩
    else:
        raise ValueError("Direction must be 'up' or 'down'")
    
    state = single_spin
    for _ in range(N - 1): 			# tensor product of N spins
        state = np.kron(state, single_spin)	# Kronecker product
        
    return state



#========================================================================
#==========================Initial state=================================

def initial_state(Ns, operators, initial_rotation):
    """
    Crea un estado inicial rotado para un sistema de N espines.
    
    Args:
        N (int): Número de espines.
        operators (list): Lista de operadores de espín (matrices numpy).
        initial_rotation (float): Ángulo de rotación inicial.
    
    Returns:
        numpy.ndarray: Vector de estado rotado.
    """

    total_operator = np.sum(operators, axis=0)			#sum operators s1 + s2 + .... 
    rotation_op = expm(-1j * initial_rotation * total_operator)	# unitary rotation
    
    psi0 = szbasis(Ns,'up')     				# initial state |0⟩^N tensor product of |0⟩ for spín)
    rotated_state = rotation_op @ psi0				# appliying rotation
    
    return rotated_state

#============================================================================
#=======================time dependent state expanded =======================
def psit_expandido(U_evals,U_evecs,coefi,t,n):
    psitt = 0
    for l in range(len(U_evals)):
        phase = U_evals[l]
        psitt += phase**n * coefi[l] * U_evecs[l]
        
    return psitt
    
#============================================================================
#=======================Fisher information =======================
def fisher(Nl, deltaBl, psitl, psit_bpl, psit_bnl):
    
    psi_deriv = (psit_bpl - psit_bnl)/(2*deltaBl)
    
    F = 4 * ( np.vdot(psi_deriv, psi_deriv) - abs(np.vdot(psitl, psi_deriv))**2)
    
    return F
    
    
#============================================================================
#=======================Function definition =======================

def f_time(tl, betal, omegal):
    print("tiempo = ", tl)
    print("beta = ", betal)
    print("omegal = ", omegal)
    return betal * np.sin(omegal * tl)
    



delta 	= 0.5
J 	= 4
N 	= 4
terms	= 2
T	= 1

theta0 	= np.pi/8
theta  	= np.pi
omega  	= np.pi/(8*T)
beta   	= 0.0
delta_b = 0.00001
dt	= 0.001
steps   = int(1e4)

np.random.seed(60)    
randomx = np.random.uniform(-delta/2, delta/2, N)
randomz = np.random.uniform(0, delta, N)
randoms = [randomx,randomz]

Hf = Hamiltonian(f_time,J,N,randoms,terms)

funcion = f_time(T, beta+delta_b, omega)
print("funccion", funcion)
print("===================================================================")

h0 = Hf(T,beta,omega)
print("h0 = ",h0)
h1 = Hf(T,0.1,omega)
print("h1 = ",h1)
h2 = Hf(T,0.5,omega)
print("h2 = ", h2)

sx = create_spin_operators(N, 'x')
psi0 = initial_state(N, sx, theta0)


def fisher_simulation(hamiltonianl, psi0l, Jl, Nl, Tl, betal, deltabl, dtl, omegal, thetal, H_termsl):
    C     = []
    sxl = create_spin_operators(Nl,'x')
    szl = create_spin_operators(Nl,'z')
    UF = build_U2(hamiltonianl, Jl, Nl, H_termsl, Tl, dtl, betal, omegal, thetal, sxl)
    UF_evals, UF_evecs = eig(UF)
    coeffs0 = [np.vdot(u_j, psi0) for u_j in UF_evecs]
    #print(UF)	
    print("coeficientes 0",coeffs0)
    UF_p1 = build_U2(hamiltonianl, Jl, Nl, H_termsl, Tl, dtl, betal - deltabl, omegal, thetal, sx)
    UFp_evals,UFp_evecs = eig(UF_p1)
    coeffsp = [np.vdot(u_j, psi0) for u_j in UFp_evecs]
    print()
    print("coeffsp",coeffsp)
    UF_m1 = build_U2(hamiltonianl, Jl, Nl, H_termsl, Tl, dtl, betal + deltabl, omegal, deltabl, sxl)
    UFm_evals,UFm_evecs = eig(UF_m1)
    coeffsm = [np.vdot(u_j, psi0) for u_j in UFm_evecs]
    print()
    print("coeffsm",coeffsm)    
    M = []
    
    fout = open("result.dat", "w")
    fout.write("# t\tmagnetization\tinfo_fisher\tmetodo_momentos\n")
    
    for n in range(1,steps):
        psit	= psit_expandido(UF_evals,UF_evecs,coeffs0,Tl,n)
        psit_p1 = psit_expandido(UFp_evals,UFp_evecs,coeffsp,Tl,n)
        psit_m1 = psit_expandido(UFm_evals,UFm_evecs,coeffsm,Tl,n)
        
        magnetization = np.vdot(psit, szl[0] @ psit)
        
        Fish = fisher(Nl, deltabl, psit, psit_p1, psit_m1)
        C.append(Fish/((2*n/np.pi)**2 * Nl) )
        M.append( magnetization)
        fout.write(f"{n:.4f}\t{magnetization.real:.16f}\t{Fish.real:.16f}\n")
        #linea_psit = " ".join(f"{z.real:.8f}+{z.imag:.8f}j" for z in psit)
        #fout.write(f"{n:.0f}\t{linea_psit}\n")
        #linea_psitp1 = " ".join(f"{z.real:.8f}+{z.imag:.8f}j" for z in psit_p1)
        #fout.write(f"{n:.0f}\t{linea_psitp1}\n")        
        #linea_psitm1 = " ".join(f"{z.real:.8f}+{z.imag:.8f}j" for z in psit_m1)
        #fout.write(f"{n:.0f}\t{linea_psitm1}\n")            
    fout.close()
    
    return C


fisher_beta2  = fisher_simulation(Hf, psi0, J, N, T, beta, delta_b, dt, omega, theta, terms)

#print(fisher_beta2)

#print()
#print(psi0)
    
    
    
    
