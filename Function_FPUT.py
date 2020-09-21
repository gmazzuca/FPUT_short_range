# Function the evolution of the FPUT chain with short range interactions

import numpy as np
import scipy as sc
import sys
import time
import matplotlib.pyplot as plt


################### MISCELLANEA ####################################################

def DHT(x):
    ''' Discrete Hartley Transform'''
    fx = np.fft.fft(x)
    return np.real(fx) - np.imag(fx)


def DHTn(x):
    ''' Discrete Hartley Transform NORMALIZED'''
    fx = np.fft.fft(x)
    return (np.real(fx) - np.imag(fx))/np.sqrt(len(x))


def log_space(t,logstep):
    ''' logarithmic scale of time '''
    
    xfine = np.arange(0,int(min([t,20])))
    while (xfine[-1] < t):
        tmp = xfine[-1]*logstep
        if tmp < t:
            xfine = np.append(xfine, tmp)
        else:
            xfine = np.append(xfine,t)

    return xfine


def periodic_difference(x,d):
    ''' give the periodic distance at d level '''

    n = len(x)
    ytmp = np.append(x[-d:],np.append(x,x[:d]))
    y = ytmp[d:] - ytmp[:-d]
    return y[d:]

def shift(A,d):
    ''' shifting elements of A's row by a their number '''
    B = np.zeros(A.shape)
    for l in range(d):
        B[l,:] = np.append(A[l,-l-1:], A[l,:-l-1])

    return B


############# Check functions (energy) #############


def periodic_energy(p,q,d,coef,chi,gamma):

    '''energy periodic toda in (p,q) variables'''
    A = np.zeros((d,len(q)))
    Acoef  = np.zeros((d,len(q)))
    for l in range(d):
        A[l,:] =  periodic_difference(q,l+1)
        Acoef[l,:] = coef[l]*A[l,:]
    qpart = np.sum( A*Acoef*(0.5 + A*(chi/3 + A*gamma*0.25)))
    ppart = 0.5*np.sum(p*p)
    
    return qpart + ppart


######### Evolution with a LP algorithm ###############
def vecp(q,d,coef,chi,gamma):
    
    ''' force for the evolution, must give dH/dq '''
    A = np.zeros((d,len(q)))
    Acoef  = np.zeros((d,len(q)))
    for l in range(d):
        A[l,:] =  periodic_difference(q,l+1)
        Acoef[l,:] = coef[l]*A[l,:]
    forcematrix =  Acoef*(1 + A*(chi + A*gamma))
    force = np.sum(- forcematrix + shift(forcematrix,d) ,0)
    return force

def leap_frog(p,q,dt,d,coef,chi,gamma):
    ''' Leap frog generic '''
    
    q_tmp = q + 0.5*dt*p
    p_new = p -  dt*vecp(q_tmp,d,coef,chi,gamma)
    q_new = q_tmp + 0.5*dt*p_new

    return(p_new,q_new)

def yo4(p,q,dt,d,coef,chi,gamma):
    ''' time step integration via yoshida4'''

    x1= 1.351207191959657
    x0 = -1.702414383919315
    
    (p1, q1) = leap_frog(p,q, x1*dt,d,coef,chi,gamma)
    (p2, q2) = leap_frog(p1,q1, x0*dt,d,coef,chi,gamma)
    (p3, q3) = leap_frog(p2,q2, x1*dt,d,coef,chi,gamma)

    return (p3,q3)

def evolution_yo4(p,q,tau,dt,d,coef,chi,gamma):
    ''' periodic evolution till time tau of the data (p,q) with a timestep  dt'''
    time = 0
    while(time < tau):
        (p,q) = yo4(p,q,dt,d,coef,chi,gamma)
        time = time + dt
        
    return(p,q)

def complete_sol(p,q,time,dt,d,coef,chi,gamma):
    ''' complete solution for the periodic case '''
    tsteps = len(time)
    particles = len(p)

    solp = np.zeros((tsteps, particles))
    solq = np.zeros((tsteps, particles))
    evostep = np.ediff1d(time)

    solp[0,:] = p
    solq[0,:] = q
    
    for k in range(tsteps - 1 ):
        (solp[k+1],solq[k+1]) = evolution_yo4(solp[k],solq[k],evostep[k], dt,d,coef, chi,gamma)

    return (solp,solq)


def circmatrix(a):
    ''' generate a circulant matrix starting from vector '''
    n = len(a)

    M = np.zeros((n,n))
    for k in range(n):
        M[k,:] = np.append(a[k:], a[:k])

    return M

def circ_root(eig):

    ''' square root of circulant matrix '''
    for k in range(len(eig)):
        if eig[k] < 0 :
            print('impossibile square root\n')
            exit()


    sqrteig = np.sqrt(eig)
    vector = DHT(sqrteig)/len(eig)

    return circmatrix(vector)



######### Starting Point - Stats #####################


def eigenvalues_force_matrix(coef,n):
    ''' eigenvalues of the interacting matrix of q '''

    d = len(coef)
    compl_coef = np.append(np.append(np.append(2*sum(coef), -coef) , np.zeros(n - 2*d-1)) , - coef[::-1])
    eig = DHT(compl_coef)
    eig[0] = 0
    return eig


def initial_condition(n,beta,eig):
    ''' Random initial conditions according the Gibbs ensemble of the UNPERTURBED chain '''
    sigma = 1/np.sqrt(beta)

    # INITIAL CONDITION ON P
    tmp_p = np.random.normal(loc=0.0, scale=sigma, size= n-1) # normal on the independent variable
    tmp_p = np.append(0,tmp_p)
    p = DHTn(tmp_p)

    # INITIAL CONDITION ON Q

    tmp_q = np.random.normal(loc=0.0, scale=sigma/np.sqrt(eig[1:]), size= n-1) # normal on the independent variable
    tmp_q = np.append(0,tmp_q)
    q = DHTn(tmp_q)

    return (p,q)

    
