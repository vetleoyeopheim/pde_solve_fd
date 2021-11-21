
"""
@author: Vetle Øye Opheim
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
    
    def burgers_eq_2d(self, init_cond, ranges, deltas, bval, h,):
        """
        Uses explicit forward steps
        Pass v=0 for inviscid fluid
        bval is a list of two lists. bval[0] is the list of boundary values for x and bval[1] for y
        """
        
        nt= int(ranges[0] / deltas[0])
        nx = int(ranges[1] / deltas[1])
        ny = int(ranges[2] / deltas[2])
        r1 = deltas[0] / deltas[1]
        r2 = deltas[0] / deltas[2]
        r3 = (deltas[0] * h) / (deltas[1] * deltas[1])
        r4 = (deltas[0] * h) / (deltas[2] * deltas[2])
        sol_x = np.zeros((nt,nx,ny))
        sol_y = np.zeros((nt,nx,ny))
        sol_x[0] = init_cond[0]
        sol_y[0] = init_cond[1]
        
        for t in range(1,nt,1):
            arr_u = np.zeros((nx,ny))
            arr_v = np.zeros((nx,ny))
            for i in range(0,nx,1):
                colmn_u = np.zeros(ny)
                colmn_v = np.zeros(ny)
                for j in range(0,ny,1):
                    if i == 0:
                        u = bval[0][0]
                        v = bval[1][0]
                    elif i == nx - 1:
                        u = bval[0][1]
                        v = bval[1][1]
                    elif j == 0:
                        u = bval[0][0]
                        v = bval[1][0]
                    elif j == ny - 1:
                        u = bval[0][1]
                        v = bval[1][1]
                    else:
                        u = sol_x[t - 1][i][j] - r1 * sol_x[t - 1][i][j] * (sol_x[t - 1][i][j] - sol_x[t - 1][i - 1][j])\
                                - r1 * sol_y[t - 1][i][j] * (sol_x[t - 1][i][j] - sol_x[t - 1][i][j - 1])\
                                + r3 * (sol_x[t - 1][i + 1][j] + sol_x[t - 1][i - 1][j] - 2 * sol_x[t - 1][i][j])\
                                + r4 * (sol_x[t - 1][i][j + 1] + sol_x[t - 1][i][j - 1] - 2 * sol_x[t - 1][i][j])
                                
                        v = sol_y[t - 1][i][j] - r2 * sol_y[t - 1][i][j] * (sol_y[t - 1][i][j] - sol_y[t - 1][i][j - 1])\
                                - r2 * sol_x[t - 1][i][j] * (sol_y[t - 1][i][j] - sol_y[t - 1][i - 1][j])\
                                + r3 * (sol_y[t - 1][i + 1][j] + sol_y[t - 1][i - 1][j] - 2 * sol_y[t - 1][i][j])\
                                + r4 * (sol_y[t - 1][i][j + 1] + sol_y[t - 1][i][j - 1] - 2 * sol_y[t - 1][i][j])
                                
                    colmn_u[j] = u
                    colmn_v[j] = v
                arr_u[i] = colmn_u
                arr_v[i] = colmn_v
            sol_x[t] = arr_u
            sol_y[t] = arr_v
        return sol_x, sol_y
