# Time_Crystals1
It is a research project coding in python and other languages like fortran90 and C++. The idea is have several version codes of simulations







La estructura de nuestro programa esta de la sigueiente forma:
A. folder Static_functions/
	statics_functions.py:
		1. sigma_i(direction)
		2. create_spin_operators(N, direction)
		3. create_spin_xyz_operators(N)
					4. disorder_x(N, delta)
					5. disorder_y(N, delta)
					6. disorder_z(N, delta)
					7. initial_state_up(N)
					8. initial_state_down(N)
					9. spin_rotator(operator, angle)
					10. check_rotator(rotator)
B. Folder hamiltonian_models/ 
    hamiltonian.py:
	    1. H_ising(J, operators)
	    2. ......
