#! /usr/bin/env python
import numpy
import pylab

"""
The Townes uses the solution to the following ODE

R'' + (1/r)*R' - R + R**3 = 0

Here we let R' = V then solve the following system using RK4 and the 
shooting method with:
V' = R - V/r - R**3
R' = V
"""

def f_V(h, r, V, R):
    if r != 0.:
        return  -V/r + R - R**3
    # Handle div by zero
    else:
        return R - R**3

def f_R(h, r, V, R):
    return V 

def k1(h, r, V, R):
    k1V = f_V(h, r, V, R)
    k1R = f_R(h, r, V, R)
    return k1V, k1R

def k2(h, r, V, R, k1V, k1R):
    k2V = f_V(h, r + .5*h, V + h*.5*k1V, R + h*.5*k1R)
    k2R = f_R(h, r + .5*h, V + h*.5*k1V, R + h*.5*k1R)
    return k2V, k2R

def k3(h, r, V, R, k2V, k2R):
    k3V = f_V(h, r + .5*h, V + h*.5*k2V, R + h*.5*k2R)
    k3R = f_R(h, r + .5*h, V + h*.5*k2V, R + h*.5*k2R)
    return k3V, k3R

def k4(h, r, V, R, k3V, k3R):
    k4V = f_V(h, r + h, V + h*k3V, R + h*k3R)
    k4R = f_R(h, r + h, V + h*k3V, R + h*k3R)
    return k4V, k4R

def step(h, r, V, R):
    k1V, k1R = k1(h, r, V, R)
    k2V, k2R = k2(h, r, V, R, k1V, k1R)
    k3V, k3R = k3(h, r, V, R, k2V, k2R)
    k4V, k4R = k4(h, r, V, R, k3V, k3R)
    
    new_V = V + (h/6.)*(k1V + 2*k2V + 2*k3V + k4V)
    new_R = R + (h/6.)*(k1R + 2*k2R + 2*k3R + k4R)

    return new_V, new_R

def RK4(h, N, IC):
    """
    Implement RK4 with N number of steps of size h with the initial 
    condition IC.
    """
    r = numpy.arange(0, N*h, h)
    R = numpy.zeros_like(r)
    V = numpy.zeros_like(r)
    
    R[0] = IC
    V[0] = 0

    for i, r_val in enumerate(r[:-1]):
        V[i+1], R[i+1] = step(h, r_val, V[i], R[i])

    return r, R

def fit(peak_x, peak_y):
    pass

if __name__=='__main__':
    x, y = RK4(0.01, 100*10, 2)
    pylab.plot(x, y)
    pylab.show()
