# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:58:02 2021

@author: Vetle
"""

class HeatEquation():
    
    def __init__(self):
        pass

    def heat_eq_3d(self, temp_0, ranges, deltas, b_vals, alpha):
        """
        Parameters:
            temp_0: Initial temperature distribution, must be 3D array
            ranges: Range of values for t,x,y,z
            deltas: deltas = [dt,dx,dy,dz]
            b_vals: boundary values [x,y,z]
            alpha: coefficient of diffusion
        """
        nt = int(ranges[0] / deltas[0])
        nx = int(ranges[1] / deltas[1])
        ny = int(ranges[2] / deltas[2])
        nz = int(ranges[3] / deltas[3])
        r1 = (alpha * deltas[0]) / (deltas[1] * deltas[1])
        r2 = (alpha * deltas[0]) / (deltas[2] * deltas[2])
        r3 = (alpha * deltas[0]) / (deltas[3] * deltas[3])
        
        temps = np.zeros((nt, nx, ny, nz))
        temps[0] = temp_0
        
        t = 1

        while t < nt - 1:
            xyz_arr = np.zeros((nx,ny,nz))
            x = 0
            while x < nx - 1:
                xy_arr = np.zeros((nx,ny))
                y = 0
                while y < ny - 1:
                    z_row = np.zeros(nz)
                    z = 0
                    while z < nz - 1:
                        if z == 0 or z == nz - 1:
                            T = b_vals[2]
                        elif y == 0 or y == ny - 1:
                            T = b_vals[1]
                        elif x == 0 or x == nx - 1:
                            T = b_vals[0]
                        else:
                            T = temps[t - 1][x][y][z] \
                            + r1 * (temps[t - 1][x - 1][y][z] + temps[t - 1][x + 1][y][z] - 2 * temps[t - 1][x][y][z]) \
                            + r2 * (temps[t - 1][x][y - 1][z] + temps[t - 1][x][y + 1][z] - 2 * temps[t - 1][x][y][z]) \
                            + r3 * (temps[t - 1][x][y][z - 1] + temps[t - 1][x][y][z + 1] - 2 * temps[t - 1][x][y][z])
                        z_row[z] = T
                        z += 1
                    xy_arr[y] = z_row
                    y += 1
                xyz_arr[x] = xy_arr
                x += 1
            temps[t] = xyz_arr
            t += 1

        return temps
    

    def heat_eq_2d(self, temp_0, ranges, deltas, b_vals, alpha):
        """
        Parameters:
            temp_0: Initial temperature distribution, must be 2D array
            ranges: Range of values for t,x,y
            deltas: deltas = [dt,dx,dy]
            b_vals: boundary values [x,y]
            alpha: coefficient of diffusion
        """
        nt = int(ranges[0] / deltas[0])
        nx = int(ranges[1] / deltas[1])
        ny = int(ranges[2] / deltas[2])
        r1 = (alpha * deltas[0]) / (deltas[1] * deltas[1])
        r2 = (alpha * deltas[0]) / (deltas[2] * deltas[2])
        
        temps = np.zeros((nt, nx, ny))
        temps[0] = temp_0
        
        t = 1
        
        while t < nt - 1:
            xy_arr = np.zeros((nx,ny))
            x = 0
            while x < nx - 1:
                y_row = np.zeros(ny)
                y = 0
                while y < ny - 1:
                    if y == 0 or y == ny:
                        T = b_vals[1]
                    elif x == 0 or x == nx:
                        T = b_vals[0]
                    else:
                        T = temps[t - 1][x][y] \
                        + r1 * (temps[t - 1][x - 1][y] + temps[t - 1][x + 1][y] - 2 * temps[t - 1][x][y]) \
                        + r2 * (temps[t - 1][x][y - 1] + temps[t - 1][x][y + 1] - 2 * temps[t - 1][x][y])
                    y_row[y] = T
                    y += 1
                xy_arr[x] = y_row
                x += 1
            temps[t] = xy_arr
            t += 1
        return temps
    
    
    def heat_eq_1d(self, temps_0, x_range, t_range, dt, dx, bval_x):
        
        """
        Parameters:
            temp_0: 1d array of initial temperatures
            x_range: range of x
            t_range: range of t
            dt: time step
            dx: space step in the x direction
            bval_x: boundary value of x
        """
        
        nx = int(x_range / dx)
        nt = int(t_range / dt)
        alpha = 0.2
        r = (alpha * dt) / (dx * dx)
        temps = np.zeros((nt,nx))
        temps[0] = temps_0
    
        
        t = 1
        
        while t < nt:
            row = np.zeros(nx)
            x = 0
            while x < nx:
                if x == 0 or x == nx - 1 :
                    temp = bval_x
                else:
                    temp = temps[t - 1][x] + r * (temps[t - 1][x - 1] + temps[t - 1][x + 1] - 2 * temps[t - 1][x])
                row[x] = temp
                x += 1
            temps[t] = row
            t += 1
        return temps
