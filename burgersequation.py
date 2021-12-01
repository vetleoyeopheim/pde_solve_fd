
"""
@author: Vetle Øye Opheim
"""

import numpy as np
from numba import njit

@njit
def burgers_eq_1d(init_cond, t_range, x_range, dx, dt, bval, v):
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


@njit
def burgers_eq_2d(init_x, init_y, ranges, deltas, bval_x, bval_y, h,):
    """
    Uses explicit forward steps
    Pass h=0 for inviscid fluid
    bval is a list of two lists. bval[0] is the list of boundary values for x and bval[1] for y
    """
    
    nt= int(ranges[0] / deltas[0])
    nx = int(ranges[1] / deltas[1])
    ny = int(ranges[2] / deltas[2])

    r1 = deltas[0] / deltas[1]
    r2 = deltas[0] / deltas[2]
    r3 = (deltas[0] * h) / (deltas[1] * deltas[1])
    r4 = (deltas[0] * h) / (deltas[2] * deltas[2])
    ux = np.zeros((nt,nx,ny))
    uy = np.zeros((nt,nx,ny))
    ux[0] = init_x[0]
    uy[0] = init_y[1]
    
    for t in range(1,nt,1):
        ux_xy = np.zeros((nx,ny))
        uy_xy = np.zeros((nx,ny))
        for i in range(0,nx,1):
            ux_y = np.zeros(ny)
            uy_y = np.zeros(ny)
            for j in range(0,ny,1):
                if i == 0:
                    u = bval_x[0]
                    v = bval_y[0]
                elif i == nx - 1:
                    u = bval_x[1]
                    v = bval_y[1]
                elif j == 0:
                    u = bval_x[2]
                    v = bval_y[2]
                elif j == ny - 1:
                    u = bval_x[3]
                    v = bval_y[3]
                else:
                    x = ux[t - 1][i][j]\
                            - r1 * ux[t - 1][i][j] * (ux[t - 1][i][j] - ux[t - 1][i - 1][j])\
                            - r1 * uy[t - 1][i][j] * (ux[t - 1][i][j] - ux[t - 1][i][j - 1])\
                            + r3 * (ux[t - 1][i + 1][j] + ux[t - 1][i - 1][j] - 2 * ux[t - 1][i][j])\
                            + r4 * (ux[t - 1][i][j + 1] + ux[t - 1][i][j - 1] - 2 * ux[t - 1][i][j])
                            
                    y = uy[t - 1][i][j]\
                            - r2 * uy[t - 1][i][j] * (uy[t - 1][i][j] - uy[t - 1][i][j - 1])\
                            - r2 * ux[t - 1][i][j] * (uy[t - 1][i][j] - uy[t - 1][i - 1][j])\
                            + r3 * (uy[t - 1][i + 1][j] + uy[t - 1][i - 1][j] - 2 * uy[t - 1][i][j])\
                            + r4 * (uy[t - 1][i][j + 1] + uy[t - 1][i][j - 1] - 2 * uy[t - 1][i][j])
                            
                ux_y[j] = x
                uy_y[j] = y
            ux_xy[i] = ux_y
            uy_xy[i] = uy_y
        ux[t] = ux_xy
        uy[t] = uy_xy
    return ux, ux
