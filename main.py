import numpy as np

from config import params
from models.pauli_matrix import build_pmatrix, sigma_i
from models.initial_state import build_initial_state
#from models.hamiltonian import build_hamiltonian
#from evolution.evolve import evolve_state
#from evolution.u_operator import operator_u
#from observables.magnetization import compute_magnetization
#from observables.fisher_information import fisher_quantum, fisher_classical
#from estimation.parameter_estimation import estimate_optimal_parameter

sx = build_pmatrix(2,'x')
print(sx[0])
print()
sx0 = sigma_i('z')
print(sx0)
print()
psi0 = build_initial_state(2,'up',sx,np.pi/8)
print(psi0)

