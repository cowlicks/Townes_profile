#! /usr/bin/env python
import numpy

"""
The Townes uses the solution to the following ODE

R'' + (1/r)*R' - R + R**3 = 0

Here we let R' = V then solve the following system using eulers method:
V' = R - V/r - R**3
R' = V
"""

def V_step(r, h, R_old, V_old):
    return V_old + h*(R_old - V_old/r - R_old**3)

def R_step(r, h, R_old, V_old):
    return R_old + h*V_old

def solver(h, N, Ri=3.0):
    """
    Solve the ODE using Euler's method.
    h is distance, Ri is R @ r=0, N is number of steps. 
    """
    r   = numpy.arange(0., h*N, h) 
    R   = numpy.array([0.]*N)
    dR  = numpy.array(R)

    R[0]    = Ri
    dR[0]   = 0
    for i in range(1, N):
        R[i]    = R[i-1] + h*dR[i-1]
        dR[i]   = dR[i-1]  + h*(R[i-1] - dR[i-1]/r[i] - R[i-1]**3)

    return R

if __name__=='__main__':
    pass
