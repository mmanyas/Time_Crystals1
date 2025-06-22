import numpy as np

from config import params
from static_functions.statics_functions import*
from hamiltonian_models.hamiltonian import* 
from evolution_time.dynamics_functions import* 
from simulations.no_averaged_simulation import*
#from observables.magnetization import compute_magnetization
#from observables.fisher_information import fisher_quantum, fisher_classical
#from estimation.parameter_estimation import estimate_optimal_parameter



N         = 2
T         = 1.0
dt        = 0.001
beta      = 0.0
delta     = 0.5
omega     = np.pi/T
theta     = np.pi
theta0    = np.pi/8
Tsteps    = int(1e4)
delta_b   = 0.0001
N_spins   = [2,4,6,8,10]
theta_eps = 0.05



simulation_no_avg(omega,theta0,theta,theta_eps,beta,delta_b,delta,T,dt,N_spins,Tsteps)

