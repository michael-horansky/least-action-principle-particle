# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:46:55 2021

@author: micha
"""

import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

m_max = 3
n_max = 3

a = 1.0
k = 25.0
        
b = a / 2.0
def monopole_v_f(x, y):
        
    #dynamic properties
    z = x / b
    s = y / b
        
    #Bound handles
    if z < 0.0:
        z = 0.0 - z
    if z >= 2.0:
        while z >= 2.0:
            z = z - 2.0
    if z > 1.0 and z < 2.0:
        z = 2.0 - z
    
    if s < 0.0:
        s = 0.0 - s
    if s >= 2.0:
        while s >= 2.0:
            s = s - 2.0
    if s > 1.0 and s < 2.0:
        s = 2.0 - s
    
    if s == 1.0 and z == 1.0:
        return(1e3)
    
    result_sum = 0.0
    for m in range(m_max):
        for n in range(n_max):
            term_1 = 1.0 / np.sqrt( (2 * m + 1 + s) ** 2 + (2 * n + 1 + z) ** 2 )
            term_2 = 1.0 / np.sqrt( (2 * m + 1 + s) ** 2 + (2 * n + 1 - z) ** 2 )
            term_3 = 1.0 / np.sqrt( (2 * m + 1 - s) ** 2 + (2 * n + 1 + z) ** 2 )
            term_4 = 1.0 / np.sqrt( (2 * m + 1 - s) ** 2 + (2 * n + 1 - z) ** 2 )
            result_sum += term_1 + term_2 + term_3 + term_4
        
    return (k / (2.0 * a) * result_sum)


my_x_min = -0.6
my_x_max = 1.6
my_y_min = -0.8
my_y_max = 0.8

x_spaces = np.linspace(my_x_min, my_x_max, 100)
y_spaces = np.linspace(my_y_min, my_y_max, 100)

v_matrix = []
for i_y in range(len(y_spaces)):
    v_matrix.append([])
    for i_x in range(len(x_spaces)):
        v_matrix[-1].append(monopole_v_f(x_spaces[i_x], y_spaces[i_y]))

plt.imshow(v_matrix, extent=[my_x_min,my_x_max,my_y_min,my_y_max])
plt.scatter([0.0, 1.0], [0.0, 0.0], s = 50, label = 'Boundary conditions')
plt.colorbar(orientation='horizontal', label = 'V [J.kg$^{-1}$]')
plt.legend()
plt.show()

#fig = plt.figure()
#ax = plt.axes(projection='3d')


#x_mesh, y_mesh = np.meshgrid(x_spaces, y_spaces)
#potential_mesh = monopole_v_f(x_mesh, y_mesh)



