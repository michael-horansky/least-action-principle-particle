# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 01:58:10 2021

@author: micha
"""


import matplotlib.pyplot as plt
import numpy as np


test_no = 10

filename = "evols/sc_1_spaces_" + str(test_no) + ".txt"

f_input = open(filename, 'r')
f_lines = f_input.readlines()
f_input.close()

evol_matrix = []

time_space = f_lines[0].split(' ')[:-1]
time_space = np.float_(time_space)
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
        

plt.title("Time evolution of trajectories of particles in scenario no. " + str(test_no))
plt.xlabel("time [s]")
plt.ylabel("displacement [m]")

for i in range(len(x_spaces)):
    if i == 0:
        my_label = 'Trivial trajectory'
    else:
        my_label = "Trajectory " + str(i)
    plt.plot(time_space, x_spaces[i], label = my_label + " x")
    plt.plot(time_space, y_spaces[i], label = my_label + " y")

plt.legend(bbox_to_anchor=(1.05, 1))
plt.show()


plt.title("Spatial trajectories of the particles in scenario no. " + str(test_no))
plt.xlabel("x [m]")
plt.ylabel("y [m]")


for i in range(len(x_spaces)):
    if i == 0:
        my_label = 'Trivial trajectory'
    else:
        my_label = "Trajectory " + str(i)
    plt.plot(x_spaces[i], y_spaces[i], label = my_label)
#plt.plot(neighbor_spaces[0], neighbor_spaces[1], label = 'Perturbed path')

plt.legend(bbox_to_anchor=(1.05, 1))
plt.show()


#  Plot the speed of each particle as a function of time

speed_spaces = []

for particle_i in range(len(x_spaces)):
    speed_spaces.append([])
    for i in range(len(time_space) - 1):
        delta_t = time_space[i + 1] - time_space[i]
        delta_x = x_spaces[particle_i][i + 1] - x_spaces[particle_i][i]
        delta_y = y_spaces[particle_i][i + 1] - y_spaces[particle_i][i]
        delta_l = np.sqrt(delta_x ** 2 + delta_y ** 2)
        speed_spaces[-1].append(delta_l / delta_t)

plt.title("Speed of the particles in scenario no. " + str(test_no))
plt.xlabel("t [s]")
plt.ylabel("|v| [m.s$^{-1}$]")


for i in range(len(x_spaces)):
    if i == 0:
        my_label = 'Trivial trajectory'
    else:
        my_label = "Trajectory " + str(i)
    plt.plot(time_space[:-1], speed_spaces[i], label = my_label)
#plt.plot(neighbor_spaces[0], neighbor_spaces[1], label = 'Perturbed path')

plt.legend(bbox_to_anchor=(1.05, 1))
plt.show()



