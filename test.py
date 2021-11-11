# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 22:04:15 2021

@author: Vetle
"""


import numpy as np
import solver

x_range = 10
t_range = 10
dx = 0.1
dt = 0.005

temps_0 = np.random.randint(0,5,(int(x_range / dx)))
temps = solver.PDESolver.heat_equation1d(temps_0, x_range, t_range, dt, dx, 1)