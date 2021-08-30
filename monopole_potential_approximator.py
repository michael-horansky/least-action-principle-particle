# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 10:54:03 2021

@author: micha
"""


import numpy as np
import datetime

m_max = 10
n_max = 10

a = 1.0
k = 1.0
        
b = a / 2.0
def monopole_v_f(my_q):
    x = my_q[0]
    y = my_q[1]
    
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
    
    if z == 1.0 and s == 1.0:
        #Inside a monopole
        return(1e20)
    result_sum = 0.0
    for m in range(m_max):
        for n in range(n_max):
            # Edit this to change the behaviour of the mesh grid
            term_1 = 1.0 / np.sqrt( (2 * m + 1 + s) ** 2 + (2 * n + 1 + z) ** 2 )
            term_2 = 1.0 / np.sqrt( (2 * m + 1 + s) ** 2 + (2 * n + 1 - z) ** 2 )
            term_3 = 1.0 / np.sqrt( (2 * m + 1 - s) ** 2 + (2 * n + 1 + z) ** 2 )
            term_4 = 1.0 / np.sqrt( (2 * m + 1 - s) ** 2 + (2 * n + 1 - z) ** 2 )
            result_sum += term_1 + term_2 + term_3 + term_4
    
    return (k / (2.0 * a) * result_sum)


f_output = open("monopole_potential_mesh.txt", "w")
mesh_density = 50
xspace = np.linspace(0.0, b, mesh_density)
yspace = np.linspace(0.0, b, mesh_density)

start_dt = datetime.datetime.now()

v_0 = monopole_v_f([0.0, 0.0])

for i_y in range(mesh_density):
    print("Calculating row", i_y + 1, "of", mesh_density)
    for i_x in range(mesh_density):
        cur_delta_v = monopole_v_f([xspace[i_x], yspace[i_y]]) - v_0
        f_output.write("%f " % cur_delta_v)
    f_output.write("\n")

f_output.close()

end_dt = datetime.datetime.now()
duration = end_dt - start_dt
print("Process terminated in", duration.total_seconds(), "seconds.")
