Townes_profile
==============

A simple numerical solution to the Townes profile in python using the 
Shooting Method and Runge-Kutta 4.

The Townes profile is given by the solution R(r) to the ODE:

    R'' + R'/r - R + R**3 = 0

With the boundary conditions:

    lim r -> inf, R -> 0
    R'(0) = 0
