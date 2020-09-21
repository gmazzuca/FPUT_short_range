# FPUT_short_range
Evolution of the classical Fermi-Pasta-Ulam-Tsigu chain with short range interaction,
i.e. the evolution of the following dynamical system:
H_FPUT = H_2(p,q) + chi*H_3(q) + gamma*H_4(q)
Where H_2(p,q) = \sum_{j=0}^{N-1} p_j^2/2 + \sum_{l=1}^m a_l(q_{j+l} - q_j)^2/2
H_3(q) = chi*\sum_{j=0}^{N-1}\sum_{l=1}^m a_l(q_{j+l} - q_j)^3/3
H_4(q) = gamma*\sum_{j=0}^{N-1}\sum_{l=1}^m a_l(q_{j+l} - q_j)^4/4 and a_l > 0

The initial date are sample at random from the Gibbs ensemble of just the quadratic part of the chain, i.e. they are sample according to:

du = e^{-beta*H_2(p,q)}dpdq/Z(beta)

where beta is the inverse of the temperature and Z(beta) is a norming constant.

To implement the evolution we used a symplectic Yoshida algorithm of order 4 implemented via a leap frog.
To work FPUT_evolution.py needs some data from the command line. In particular:
 - inverse_temperature: the inverse of the temperature of the system, i.e. beta;
 - final_time: final time of the evolution;
 - number_particles: number of particles of the chain;
 - chi, gamma: values for the two constants;
 - integration step, usually between 0.01 and 0.1;
 - log step: logarithmic scale for the evolution;


In Function_FPUT.py there are all the functions needed by FPUT_evolution.py to actually compute the evolutions of the chain.
