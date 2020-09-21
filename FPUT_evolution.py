# Evolution of the FPUT chain with short range interaction using a Yoshida algorothm of order 4
# The system is described  by the Hamiltonian H_FPUT = H_2 + chi*H_3 + gamma*H_4
# Where H_2 = \sum_{j=0}^{N-1} p_j^2/2 + \sum_{l=1}^m a_l(q_{j+l} - q_j)^2/2
# H_3 = chi*\sum_{j=0}^{N-1}\sum_{l=1}^m a_l(q_{j+l} - q_j)^3/3
# H_4 = gamma*\sum_{j=0}^{N-1}\sum_{l=1}^m a_l(q_{j+l} - q_j)^4/4
# a_l > 0
# we will refer to a_l as the springs strenght 

import numpy as np
import scipy as sc
import sys
import matplotlib.pyplot as plt
import Function_FPUT as fpu


np.random.seed()
if (len(sys.argv) < 7):
    print('error: give inverse_temperature final_time number_particles chi gamma integration_step log_step\n')
    exit()

# PARAMETERS FOR EVOLUTION FROM COMMAND LINE
beta = int(sys.argv[1])
final_time = float(sys.argv[2])
number_particles = int(sys.argv[3])
chi = float(sys.argv[4])
gamma = float(sys.argv[5])
integration_step = float(sys.argv[6])
log_step = float(sys.argv[7])


springs = np.array([1,0.125, 7/72]) #vector of the springs strenght
d = len(springs)


timer = fpu.log_space(final_time,log_step) # time sample
numstrobe = len(timer)
    
############ The game start ##################
eig_force_matrix = fpu.eigenvalues_force_matrix(springs,number_particles) # eigeenvalues for the matrix
M = fpu.circ_root(eig_force_matrix)
for k in range(number_particles):
    if eig_force_matrix[k] < 0 :
        print('error!! negative eigenvalue')
        exit()
# Evolution
(p,q) = fpu.initial_condition(number_particles,beta,eig_force_matrix)
(psol,qsol) = fpu.complete_sol(p,q,timer,integration_step,d,springs, chi, gamma)

print(psol,qsol)
    
    
    

