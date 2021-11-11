# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:58:02 2021

@author: Vetle
"""


import numpy as np



class PDESolver():
    
    def __init__(self):
        pass


    def heat_equation2d(temp_0, x_range, y_range, t_range, dt, dx, dy):
        """TODO """
        nx = int(x_range / dx)
        ny = int(y_range / dy)
        nt = int(t_range / dt)
        
        alpha = 0.2
        
        temps = 0
        
        return temps
    
    
    def heat_equation1d(temps_0, x_range, t_range, dt, dx, bval_x):
        
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


