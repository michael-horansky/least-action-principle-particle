# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 12:51:01 2021

@author: micha
"""

import matplotlib.pyplot as plt
import numpy as np


test_no = 10

filename = "evols/sc_1_a_g_evol_" + str(test_no) + ".txt"

f_input = open(filename, 'r')
f_lines = f_input.readlines()
f_input.close()

evol_matrix = []

for i in range(len(f_lines)):
    evol_matrix.append(f_lines[i].split(' ')[:-1])
    evol_matrix[-1] = np.float_(evol_matrix[-1])

plt.title("Evolution of $\\nabla S$ in M-H alg. - scenario no. " + str(test_no))
plt.xlabel("Iteration")
plt.ylabel("$\\nabla S$")

for i in range(len(evol_matrix)):
    plt.plot(evol_matrix[i], label = "Trajectory " + str(i + 1))

plt.legend(bbox_to_anchor=(1.05, 1))
plt.show()