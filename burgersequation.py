# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 18:56:41 2021

@author: Vetle
"""

import numpy as np
import matplotlib.pyplot as plt

class BurgersEquation:
    
    """
    Solves the Burgers equation in 1D and 2D using a finite differences method
    Contains solvers for cases with both viscosity and inviscid fluids
    1D equation: du/dt + u * du/dx = v * d2u/dx2
    """
    
    def __init__(self):
        pass
    
    def burgers_eq_1d(self, init_cond, t_range, x_range, dx, dt, bval, v):
        """Pass v=0 for inviscid fluid"""
        
        nx = int(x_range / dx)
        nt = int(t_range / dt)
        r1 = dt / dx
        r2 = (dt * v) / (dx * dx)
        sol = np.zeros((nt,nx))
        sol[0] = init_cond
        
        for t in range(1,nt,1):
            row = np.zeros(nx)
            for x in range(0,nx,1):
                if x == 0:
                    u = bval[0]
                elif x == nx - 1 :
                    u = bval[1]
                else:
                    u = sol[t - 1][x] - r1 * sol[t - 1][x] * (sol[t - 1][x] - sol[t - 1][x - 1]) + \
                        r2 * (sol[t - 1][x + 1] + sol[t - 1][x - 1] - 2 * sol[t - 1][x])
                row[x] = u
            sol[t] = row
        return sol
    
    def burgers_eq_2d(self, init_cond, ranges, deltas, bval, v):
        """
        Pass v=0 for inviscid fluid
        bval is a list of two lists. bval[0] is the list of boundary values for x and bval[1] for y
        """
        
        nt = int(ranges[0] / deltas[0])
        nx = int(ranges[1] / deltas[1])
        ny = int(ranges[2] / deltas[2])
        r1 = deltas[0] / deltas[1]
        r2 = deltas[0] / deltas[2]
        r3 = (deltas[0] * v) / (deltas[1] * deltas[1])
        r4 = (deltas[0] * v) / (deltas[2] * deltas[2])
        sol = np.zeros((nt,nx,ny))
        sol[0] = init_cond
        
        for t in range(1,nt,1):
            arr = np.zeros((nx,ny))
            for x in range(0,nx,1):
                row = np.zeros(ny)
                for y in range(0,ny,1):
                    if x == 0:
                        u = bval[0][0]
                    elif y == 0:
                        u = bval[1][0]
                    elif x == nx - 1 :
                        u = bval[0][1]
                    elif y == ny - 1:
                        u = bval[1][1]
                    else:
                        u = sol[t - 1][x] - r1 * sol[t - 1][x] * (sol[t - 1][x] - sol[t - 1][x - 1]) + \
                            r2 * (sol[t - 1][x + 1] + sol[t - 1][x - 1] - 2 * sol[t - 1][x])
                    row[y] = u
                arr[x] = row
            sol[t] = row
        return sol
        
"""
dx = 0.25
dt = 0.01
t_range = 10
x_range = 10
bval = [0,1]
v = 0
init_cond = np.ones(int(x_range / dx))
be = BurgersEquation()
sol = be.burgers_eq_1d(init_cond, t_range, x_range, dx, dt, bval, v)
plt.plot(sol[950])
"""

ranges = [10,10,10]
deltas = [0.01,0.2,0.2]
v = 0
bval = [[0,1],[0,1]]
init_cond = np.ones((int(ranges[1] / deltas[1]), int(ranges[2] / deltas[2])))
    
be = BurgersEquation()
sol = be.burgers_eq_2d(init_cond, ranges, deltas, bval, v)
plt.plt(sol[100][10])