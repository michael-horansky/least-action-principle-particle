# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 16:46:00 2021

@author: micha
"""


#  This is just a fun thing to simulate a particle using its obtained
#  initial position and velocity to see how much its trajectory actually
#  differs from what the LAP algorithm says.

import numpy as np
import matplotlib.pyplot as plt

#---- potential function babyyy (i should put this in a separate file and
#---- everything should just refere to it cuz these parameters are crazy)

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

def acceleration(x, y):
    h = 0.01
    a_x = - (monopole_v_f(x + h, y) - monopole_v_f(x - h, y)) / (2.0 * h)
    a_y = - (monopole_v_f(x, y + h) - monopole_v_f(x, y - h)) / (2.0 * h)
    return(a_x, a_y)

def runge_kutta(y_n, t_n, f_y, dt):
    k_1 = f_y(t_n, y_n)
    k_2 = f_y(t_n + dt / 2.0, y_n + dt * k_1 / 2.0)
    k_3 = f_y(t_n + dt / 2.0, y_n + dt * k_2 / 2.0)
    k_4 = f_y(t_n + dt, y_n + dt * k_3)
    result = y_n + (1.0 / 6.0) * dt * (k_1 + 2 * k_2 + 2 * k_3  + k_4)
    return(result)


#---- First load the LAP trajectories cuz why tf not

test_no = 10

filename = "evols/sc_1_spaces_" + str(test_no) + ".txt"

f_input = open(filename, 'r')
f_lines = f_input.readlines()
f_input.close()

evol_matrix = []

time_space = f_lines[0].split(' ')[:-1]
time_space = np.float_(time_space)
v_x_0_list = []
v_y_0_list = []
x_spaces = []
y_spaces = []

for i in range(1, len(f_lines)):
    evol_matrix.append(f_lines[i].split(' ')[:-1])
    x_spaces.append([])
    y_spaces.append([])
    for j in range(len(evol_matrix[-1])):
        evol_matrix[-1][j] = evol_matrix[-1][j].split(',')
        evol_matrix[-1][j] = np.float_(evol_matrix[-1][j])
        x_spaces[-1].append(evol_matrix[-1][j][0])
        y_spaces[-1].append(evol_matrix[-1][j][1])
    delta_t = time_space[1] - time_space[0]
    delta_x = evol_matrix[-1][1][0] - evol_matrix[-1][0][0]
    delta_y = evol_matrix[-1][1][1] - evol_matrix[-1][0][1]
    v_x_0_list.append(delta_x / delta_t)
    v_y_0_list.append(delta_y / delta_t)

"""for i in range(1, len(f_lines)):
    evol_matrix.append(f_lines[i].split(' ')[:-1])
    evol_matrix[-1][0] = evol_matrix[-1][0].split(',')
    evol_matrix[-1][0] = np.float_(evol_matrix[-1][0])
    evol_matrix[-1][1] = evol_matrix[-1][1].split(',')
    evol_matrix[-1][1] = np.float_(evol_matrix[-1][1])
    delta_t = time_space[1] - time_space[0]
    delta_x = evol_matrix[-1][1][0] - evol_matrix[-1][0][0]
    delta_y = evol_matrix[-1][1][1] - evol_matrix[-1][0][1]
    v_x_0_list.append(delta_x / delta_t)
    v_y_0_list.append(delta_y / delta_t)"""
print("A love supremum, a function infimum.")
def simulate_trajectory(x_0, y_0, v_x_0, v_y_0, t_max):
    
    dt = 0.001
    t = 0.0
    x = x_0
    y = y_0
    v_x = v_x_0
    v_y = v_y_0
    traj_x = []
    traj_y = []
    while t < 0.25:
        # runge kutta pls bro
        a_x, a_y = acceleration(x, y)
        
        new_x = x + dt * v_x + (dt * dt / 2.0) * a_x
        new_y = y + dt * v_y + (dt * dt / 2.0) * a_y
        new_v_x = v_x + dt * a_x
        new_v_y = v_y + dt * a_y
        
        x = new_x
        y = new_y
        v_x = new_v_x
        v_y = new_v_y
        
        traj_x.append(x)
        traj_y.append(y)
        
        t += dt
        #print(t < 2.0)
    
    return(traj_x, traj_y)

x_spaces_sim = []
y_spaces_sim = []
my_t = 2.0
for i in range(len(x_spaces)):
    print("Simulating trajectory no.", i, "of", len(x_spaces) - 1)
    cur_x_spaces, cur_y_spaces = simulate_trajectory(x_spaces[i][0], y_spaces[i][0], v_x_0_list[i], v_y_0_list[i], my_t)
    x_spaces_sim.append(cur_x_spaces)
    y_spaces_sim.append(cur_y_spaces)

which_traj = 5

plt.plot(x_spaces[which_traj],     y_spaces[which_traj],     label = 'Alg. result')
plt.plot(x_spaces_sim[which_traj], y_spaces_sim[which_traj], label = 'Sim. result')



plt.legend()
plt.show()

#print(v_y_0_list)

