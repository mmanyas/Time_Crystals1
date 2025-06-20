import numpy as np

from config import params
from static_functions.statics_functions import*
#from models.pauli_matrix import build_pmatrix, sigma_i
#from models.initial_state import build_initial_state
#from models.hamiltonian import build_hamiltonian
#from evolution.evolve import evolve_state
#from evolution.u_operator import operator_u
#from observables.magnetization import compute_magnetization
#from observables.fisher_information import fisher_quantum, fisher_classical
#from estimation.parameter_estimation import estimate_optimal_parameter


sx,_, sz = create_spin_xyz_operators(3)

print(disorder_x(2, 0.5))
