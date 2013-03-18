#! /usr/bin/env python
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

def solver(r, h, ri=1, N):
    """
    Solve the ODE using Euler's method.
    h is distance, ri is R @ r=0, N is number of steps
    """
    pass

if __name__=='__main__':
    pass
