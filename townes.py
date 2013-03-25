#! /usr/bin/env python
import numpy
import pylab

H   = 0.001     # RK4 step size
N   = 1000*9*2.3   # RK4 step number

"""
The Townes uses the solution to the following ODE

R'' + (1/r)*R' - R + R**3 = 0

Here we let R' = V then solve the following system using RK4 and the 
shooting method with:
V' = R - V/r - R**3
R' = V
"""

# Solving R with Runge-Kutta
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

# Finding the best R(0) with the shooting method. The root finding is 
# done with the secant method.
def shoot(IC):
    x, y    = RK4(H, N, IC)
    val     = y[-1]
    return val

def secant_method(domain, f, tol, max_iter):
    i = domain[0]
    i_old = domain[-1]
    y = f(i)
    count = 0
    while abs(y) > tol:
        i_new = i - f(i)*(i - i_old)/(f(i) - f(i_old))
        y = f(i_new)
        i_old = i
        i = i_new
        count += 1
        print "count = " + str(count)
        print "y = " + str(y)
        print "ic = " + str(i)
        if count > max_iter:
            print "max_iter exceeded"
            print counts
    return i

# Fitting the profile to data.
def fit(data, townes):
    """
    Fit the Townes profile to the given data. Such that:
    R(r) => a*R(a*r)
    """
    dx  = data[0]
    dy  = data[1]
    tx  = townes[0]
    ty  = townes[1]

    # Fit the y scale.
    dmax    = max(dy)
    tmax    = max(ty)
    ty      = ty * dmax/tmax

    # Fit the x scale.
    half_dmax_indx  = min(range(len(dy)), key=lambda i: abs(dy[i] - dmax/2))
    half_tmax_indx  = min(range(len(ty)), key=lambda i: abs(ty[i] - dmax/2))
    tx      = tx * dx[half_dmax_indx]/tx[half_tmax_indx] 
    return tx, ty

if __name__=='__main__':
    # best ic so far:
    # by hand: 2.20620158
    # shooting: 2.20620157567
    ic  = secant_method([2.2062015, 2.2062016], shoot, 0.0000001, 100)
    x, y = RK4(H, N, ic)
    pylab.plot(x, y)
    pylab.show()
