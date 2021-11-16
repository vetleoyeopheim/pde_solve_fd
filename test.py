"""
@author: Vetle
"""


import numpy as np
import heatequation

x_range = 10
t_range = 10
dx = 0.1
dt = 0.005

temps_0 = np.random.randint(0,5,(int(x_range / dx)))
temps = solver.PDESolver.heat_equation1d(temps_0, x_range, t_range, dt, dx, 1)
