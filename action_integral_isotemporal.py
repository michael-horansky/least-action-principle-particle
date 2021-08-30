# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 15:57:53 2021

@author: micha
"""

import numpy as np
from scipy.interpolate import CubicSpline, interp1d
import scipy.integrate as integrate
import scipy.optimize as optimize
import matplotlib.pyplot as plt

from pynverse import inversefunc
import datetime

#Physical quantities

T_0 = 10.0

#Debug functions

#Free fall solver
def free_fall(q_0, t_0, q_1, t_1, my_timespace):
    delta_x = q_1[0] - q_0[0]
    delta_y = q_1[1] - q_0[1]
    delta_t = t_1 - t_0
    v_x_0 = delta_x / delta_t
    v_y_0 = (delta_y + 0.5 * g * delta_t * delta_t) / delta_t
    
    x_space = np.zeros(len(my_timespace))
    y_space = np.zeros(len(my_timespace))
    
    for i in range(len(my_timespace)):
        x_space[i] = q_0[0] + v_x_0 * my_timespace[i]
        y_space[i] = q_0[1] + v_y_0 * my_timespace[i] - 0.5 * g * my_timespace[i] * my_timespace[i]
    return  x_space, y_space


#General functions
def vector_distance(vec1, vec2):
    D = 0.0
    for i in range(len(vec1)):
        D += (vec1[i] - vec2[i]) ** 2
    return(np.sqrt(D))

def vector_size(vec):
    D = 0.0
    for i in range(len(vec)):
        if type(vec[i]) == type([]):
            D += vector_size(vec[i]) ** 2
        else:
            D += vec[i] ** 2
    return(np.sqrt(D))

def vector_scaling(vec, scalar):
    result = []
    for item in vec:
        result.append(item * scalar)
    return(result)

def vector_sum(vec1, vec2):
    result_vec = []
    for i in range(len(vec1)):
        result_vec.append(vec1[i] + vec2[i])
    return(result_vec)

def vector_difference(vec1, vec2):
    return(vector_sum(vec1, vector_scaling(vec2, -1.0)))


class LAP:
    def __init__(self, my_q_0, my_t_0, my_q_1, my_t_1, my_m, my_DoF, my_og_potential):
        #Properties of the simulation
        self.x_min = 50.0
        self.x_max = -50.0

        self.delta_t = 1e-2#time step = epsilon
        self.delta_x = 1e-2

        #Conditions for the trajectory
        #stc = spacetime coordinate = [q, t]
        self.q_0 = my_q_0
        self.q_1 = my_q_1
        self.t_0 = my_t_0
        self.t_1 = my_t_1
        
        self.m = my_m
        
        self.spatial_D = len(self.q_0)
        
        self.DoF = my_DoF
        self.q_t = []
        
        self.og_potential = my_og_potential
        self.potential = self.og_potential(self)
        
        #print(type(self.potential))
        
        #Create trajectory
        #  DoF = degrees of freedom = number of tweakable points
        #  First approximation = linear betweeen q_0 and q_1
        self.trajectory = np.linspace(self.q_0, self.q_1, self.DoF + 2)
        self.timespace = np.linspace(self.t_0, self.t_1, self.DoF + 2)
        #self.trajectory_interpolate()
    
    #Copy init used in calculating the gradient and finding neighbouring trajectories
    @classmethod
    def copy(cls, other_LAP):
        
        my_copy = cls(other_LAP.q_0.copy(), other_LAP.t_0, other_LAP.q_1.copy(), other_LAP.t_1, other_LAP.m, other_LAP.DoF, other_LAP.og_potential)
        #Properties of the simulation
        my_copy.x_min = other_LAP.x_min
        my_copy.x_max = other_LAP.x_max

        my_copy.delta_t = other_LAP.delta_t#time step = epsilon
        my_copy.delta_x = other_LAP.delta_x

        #Conditions for the trajectory
        #stc = spacetime coordinate = [q, t]
        my_copy.spatial_D = other_LAP.spatial_D
        
        
        #Create trajectory
        my_copy.trajectory = other_LAP.trajectory.copy()
        my_copy.timespace = other_LAP.timespace.copy()
        #my_copy.trajectory_interpolate()
        return(my_copy)
            
    def trajectory_interpolate(self):
        self.q_t = []
        for spatial_i in range(self.spatial_D):
            self.q_t.append(interp1d(self.timespace, self.trajectory[:,spatial_i], kind="cubic"))

    def get_q(self, my_t, convert_to_number = False):
        if type(my_t) == type(np.array([])):
            if self.spatial_D == 1:
                result_space = np.zeros(len(my_t))
                for i in range(len(my_t)):
                    result_space[i] = self.get_q(my_t[i], True)
                return(result_space)
            else:
                result_space = []
                for i in range(len(my_t)):
                    result_space.append(self.get_q(my_t[i]))
                return(result_space)
        
        #assume local oddity about the endpoints
        if my_t < self.t_0:
            return(vector_scaling(self.q_0, 2.0) - self.get_q(2.0 * self.t_0 - my_t, convert_to_number))
        if my_t > self.t_1:
            return(vector_scaling(self.q_1, 2.0) - self.get_q(2.0 * self.t_1 - my_t, convert_to_number))
        
        
        cur_q = np.zeros(self.spatial_D)
        for spatial_i in range(self.spatial_D):
            cur_q[spatial_i] = self.q_t[spatial_i](my_t)
        if self.spatial_D == 1 and convert_to_number:
            return(cur_q[0])
        return(cur_q)
    
    def get_path_geometry(self, N = 100, return_timespace = False):
        #A tool for generating the geometry of a multidimensional path
        path_timespace = np.linspace(self.t_0, self.t_1, N)
        q_space = self.get_q(path_timespace)
        if return_timespace:
            return(list(zip(*q_space)), path_timespace)
        return(list(zip(*q_space)))
    
    #def potential(self, my_q):
    #    height = my_q[1]
    #    return (self.m * height * g)
    
    def get_lagrangian(self, my_t):
        speed = vector_distance(self.get_q(my_t + self.delta_t), self.get_q(my_t - self.delta_t)) / (2.0 * self.delta_t)
        T = 0.5 * self.m * speed * speed
        U = self.potential(self.get_q(my_t))
        return(T - U)
    
    def get_action(self):
        result, abserr = integrate.quad(self.get_lagrangian, self.t_0, self.t_1)
        return(result)
    
    def get_action_gradient(self):
        result_gradient = []
        for i in range(self.DoF):
            trajectory_index = i + 1
            result_gradient.append([])
            for j in range(self.spatial_D):
                tweaked_path = LAP.copy(self)
                
                delta_position = np.zeros(self.spatial_D)
                delta_position[j] = self.delta_x
                
                tweaked_path.trajectory[trajectory_index] += delta_position
                tweaked_path.trajectory_interpolate()
                action_plus = tweaked_path.get_action()
                tweaked_path.trajectory[trajectory_index] -= 2.0 * delta_position
                tweaked_path.trajectory_interpolate()
                action_minus = tweaked_path.get_action()
            
                partial_derivative = (action_plus - action_minus) / (2.0 * self.delta_x)
                result_gradient[-1].append(partial_derivative)
        return(result_gradient)
    
    def get_gradient_gradient(self):
        result_gradient = []
        for i in range(self.DoF):
            trajectory_index = i + 1
            result_gradient.append([])
            for j in range(self.spatial_D):
                tweaked_path = LAP.copy(self)
                
                delta_position = np.zeros(self.spatial_D)
                delta_position[j] = self.delta_x
                
                tweaked_path.trajectory[trajectory_index] += delta_position
                tweaked_path.trajectory_interpolate()
                action_plus = vector_size(tweaked_path.get_action_gradient())
                tweaked_path.trajectory[trajectory_index] -= 2.0 * delta_position
                tweaked_path.trajectory_interpolate()
                action_minus = vector_size(tweaked_path.get_action_gradient())
            
                partial_derivative = (action_plus - action_minus) / (2.0 * self.delta_x)
                result_gradient[-1].append(partial_derivative)
        return(result_gradient)
    
    #--------------- OPTIMIZATION ALGORITHMS ------------------------
    
    #---- Get the trivial solution using conservation of energy
    #  We know the geometry of the trivial solution - that is a straight line
    #  Hence we know v(x), and by varying v_0 we want to obtain t(x_1) = t_1
    #  Then, the perturbation will just introduce a rotation to the geometry
    #  (or, equivalently, a rotation to the potential function)
    
    #  Since t_{x_1}(v_0) is a monotone function, we can use Newton's method.
    #  The function assumes that the geometry of the correct trajectory
    #  between the endpoints is a straight line
    
    def get_trivial_solution(self, triv_q_0, triv_t_0, triv_q_1, triv_t_1, N = 50, init_guess = -1.0, method = 'scipy'):
        total_distance = vector_distance(triv_q_1, triv_q_0)
        delta_q = vector_sum(triv_q_1, vector_scaling(triv_q_0, -1.0))
        if init_guess == -1.0:
            init_guess = total_distance / (triv_t_1 - triv_t_0)
        actual_v_0 = init_guess
        #Parametrize the trajectory
        def cur_q(x):
            displacement = vector_scaling(delta_q, x / total_distance)
            return(vector_sum(triv_q_0, displacement))
        
        def cur_v(x, v_0, power = 1):
            sqrt_inside = v_0 ** 2 - 2 / self.m * (self.potential(cur_q(x)) - self.potential(triv_q_0))
            if sqrt_inside > 0:
                v = np.sqrt(sqrt_inside)
            else:
                v = 1.0e-10
            if power == 1:
                return(v)
            if power == -1:
                return(1 / v)
            if power == -2:
                return(1 / (v ** 2))
        
        
        def t_func(v_0):
            #print("Testing:", v_0)
            value, num_error = integrate.quad(lambda x: cur_v(x, v_0, -1), 0.0, total_distance)
            return(value)
        
        def t_x(my_x, v_0):
            value, num_error = integrate.quad(lambda x: cur_v(x, v_0, -1), 0.0, my_x)
            return(value)
        
        if method == 'newton':
            current_guess = init_guess
            for i in range(N):
                denominator, den_error = integrate.quad(cur_v, 0.0, total_distance, (current_guess, -2))
                new_guess = current_guess - ((triv_t_1 - triv_t_0) - t_func(current_guess)) / denominator
                #try abs_val, as speed cannot be negative
                new_guess = np.abs(new_guess)
                current_guess = new_guess
        
            #return(current_guess)
            actual_v_0 = current_guess
        
        if method == 'scipy':
            #max_V_x = optimize.fminbound(lambda x: -self.potential(cur_q(x)), 0.0, total_distance)
            #max_V = 1.001 * self.potential(cur_q(max_V_x))
            #max_V = 10.0
            #min_v = np.sqrt(2.0 / self.m * max_V)
        
            #cur_xspace = np.linspace(0.0, total_distance, 100)
            #plt.plot(cur_xspace, cur_v(cur_xspace, min_v, -1))
        
            #print(max_V)
            #print("Min v =", min_v)
            #print("Time of min v =", t_func(min_v))
            actual_v_0 = optimize.root(lambda fx: (triv_t_1 - triv_t_0) - t_func(fx), init_guess).x[0]
            #return(optimize.root(t_func, init_guess).x[0])
        actual_v_0 = np.abs(actual_v_0)
        #Now find x(t) as an inverse of t(x)
        #cur_xspace = np.linspace(0.0, total_distance, 10)
        #cur_tspace = t_x(cur_xspace, actual_v_0)
        #print(self.timespace)
        #print(self.timespace[:-1])
        #x_t = inversefunc(t_x, self.timespace, args=(actual_v_0))
        for i in range(1, self.DoF + 1):
            cur_x_t = inversefunc(t_x, self.timespace[i], args=(actual_v_0))
            self.trajectory[i] = cur_q(cur_x_t)
            #print(cur_x_t)
            
        print("Optimization ended. Current t_1 =", t_func(actual_v_0))
        return(actual_v_0)
    
    #  Once we have the trivial solution, we can now introduce a perturbation.
    #  Either change the potential function by accessing the LAP instance's
    #  properties (to perturb the potential), or use the following function
    #  to change the trajectory slightly after you initialize the trivial
    #  solution (to perturb the end conditions). This serves as a reasonable
    #  initial guess for the actual pertrubed solution.
    #  Note: you can't perturb the time conditions. That would be unreasonable.
    
    def perturb_conditions(self, new_q_0, new_q_1):
        delta_q_old = vector_difference(self.q_1, self.q_0)
        delta_q_new = vector_difference(new_q_1 , new_q_0 )
        
        #delta_t_vec = vector_scaling(vector_difference(delta_q_1, delta_q_0), 1.0 / (1.0 + self.DoF))
        for i in range(0, self.DoF + 2):
            k_n = vector_size(vector_difference(self.trajectory[i], self.q_0)) / vector_size(delta_q_old)
            self.trajectory[i] = vector_sum(new_q_0, vector_scaling(delta_q_new, k_n))
        self.q_0 = new_q_0
        self.q_1 = new_q_1
        
    #---- Gradient descent
    
    def gradient_descent(self, N=30, h=0.03):
        self.trajectory_interpolate()
        for i in range(N):
            cur_gradient = self.get_action_gradient()
            #cur_gradient = self.get_action_gradient()
            for j in range(self.DoF):
                self.trajectory[j + 1] = vector_sum(self.trajectory[j + 1], vector_scaling(cur_gradient[j], - h))
                #TODO -- what if 2D?
            self.trajectory_interpolate()
            #print(cur_gradient)
            print("Iteration", i, "with action gradient squared =", vector_size(cur_gradient))

    #---- Metropolis - Hastings algorithm
    
    def find_neighbor(self, my_std):
        neighbor = LAP.copy(self)
        for i in range(neighbor.DoF):
            for j in range(neighbor.spatial_D):
                neighbor.trajectory[i + 1][j] = np.random.normal(self.trajectory[i + 1][j], my_std)
        return(neighbor)
    
    def metropolis_hastings(self, N=30, std = 0.1, bias = 1.0):
        print("Initializing the Metropolis-Hastings algorithm:")
        self.trajectory_interpolate()
        
        cur_action_gradient = vector_size(self.get_action_gradient())
        a_g_evol = [cur_action_gradient]
        for i in range(N):
            print("  Iteration", i+1, "of", N)
            neighbor = self.find_neighbor(std)
            neighbor.trajectory_interpolate()
            neighbor_action_gradient = vector_size(neighbor.get_action_gradient())
            print("    My = %.3f, neighbor's = %.3f" % (cur_action_gradient, neighbor_action_gradient))
            if(neighbor_action_gradient < cur_action_gradient):
                acceptance_probability = 1.0
            else:
                acceptance_probability = np.exp( - bias * neighbor_action_gradient / cur_action_gradient)
            print("    A = %.3f" % acceptance_probability)
            if np.random.rand() < acceptance_probability:
                #accept
                self.trajectory = neighbor.trajectory
                self.trajectory_interpolate()
                cur_action_gradient = vector_size(self.get_action_gradient())
                #cur_action_gradient = neighbor_action_gradient
                print("    Accepted")
            a_g_evol.append(cur_action_gradient)
        return(a_g_evol)


#---- Useful procedures
def compare_DoF(my_q_0, my_t_0, my_q_1, my_t_1, my_m, starting_DoF, ending_DoF, N=30):
    compare_spaces = []
    for i in range(starting_DoF, ending_DoF + 1):
        cur_LAP = LAP(my_q_0, my_t_0, my_q_1, my_t_1, my_m, i)
        print("Checking LAP", i - starting_DoF + 1, "of", ending_DoF + 1 - starting_DoF)
        cur_LAP.gradient_descent(N)
        compare_spaces.append(cur_LAP.get_path_geometry())
    return(compare_spaces)

def plot_paths_DoF(spaces, my_plt, starting_DoF):
    for i in range(len(spaces)):
        my_plt.plot(spaces[i][0], spaces[i][1], label = str(starting_DoF + i) + ' degrees of freedom')


#my_array = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
#rot_array = list(zip(*my_array))
#print(rot_array)

#Properties of the simulation
g = 10.0
def grav_v(LAP_instance):
    def grav_v_f(my_q):
        height = my_q[1]
        return (LAP_instance.m * height * g)
    return(grav_v_f)
def sin_v(LAP_instance):
    def sin_v_f(my_q):
        height = my_q[0]
        return (4.0 * np.sin(height))
    return(sin_v_f)

def monopole_grid_v(LAP_instance):
    #constants
    m_max = 3
    n_max = 3

    a = 1.0
    k = 25.0
        
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
        
        result_sum = 0.0
        for m in range(m_max):
            for n in range(n_max):
                term_1 = 1.0 / np.sqrt( (2 * m + 1 + s) ** 2 + (2 * n + 1 + z) ** 2 )
                term_2 = 1.0 / np.sqrt( (2 * m + 1 + s) ** 2 + (2 * n + 1 - z) ** 2 )
                term_3 = 1.0 / np.sqrt( (2 * m + 1 - s) ** 2 + (2 * n + 1 + z) ** 2 )
                term_4 = 1.0 / np.sqrt( (2 * m + 1 - s) ** 2 + (2 * n + 1 - z) ** 2 )
                result_sum += term_1 + term_2 + term_3 + term_4
        
        return (k / (2.0 * a) * result_sum)
    return(monopole_v_f)

#---- Optimized potential: mesh reader ----
f_monopole_mesh = open("monopole_potential_mesh.txt", "r")
monopole_mesh_lines = f_monopole_mesh.readlines()
f_monopole_mesh.close()
monopole_mesh = []
for i_y in range(len(monopole_mesh_lines)):
    monopole_mesh.append(monopole_mesh_lines[i_y].split(" ")[:-1])
    monopole_mesh[-1] = np.float_(monopole_mesh[-1])

def interpolated_monopole_grid_v(LAP_instance):
    #constants
    
    m_max = 10
    n_max = 10

    a = 1.0
    k = 3.0
        
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
        
        #Calculate the mesh indexes
        delta_z = 1.0 / (len(monopole_mesh[0]) - 1.0)
        delta_s = 1.0 / (len(monopole_mesh) - 1.0)
        n_x = int(np.floor(z / delta_z))
        n_y = int(np.floor(s / delta_s))
        if n_x == len(monopole_mesh[0]) - 1:
            n_x_1 = n_x
        else:
            n_x_1 = n_x + 1
        if n_y == len(monopole_mesh) - 1:
            n_y_1 = n_y
        else:
            n_y_1 = n_y + 1
        
        #remainders
        r_x = z - n_x * delta_z
        r_y = s - n_y * delta_s
        #bilinear interpolation
        result_v = (1 - r_x) * (1 - r_x) * monopole_mesh[n_y][n_x] + r_x * (1 - r_y) * monopole_mesh[n_y][n_x_1] + (1 - r_x) * r_y * monopole_mesh[n_y_1][n_x] + r_x * r_y * monopole_mesh[n_y_1][n_x_1]
        #print("Result =", result_v)
        return(result_v)
    return(monopole_v_f)

#print()

h = 0.001
dx = h

my_t0 = 0.0
my_t1 = 0.25


#---- Scenario 1: Monopole potential grid, short path ----

print("---- Scenario 1: Monopole potential grid, short path ----")
sc_1_a = 1.0
sc_1_m = 1.0

sc_1_a_c = 1.0

start_dt = datetime.datetime.now()





my_LAP = LAP([0.0, 0.0], my_t0, [sc_1_a_c * sc_1_a, 0.0], my_t1, sc_1_m, 4, monopole_grid_v)
print("v_0 =", my_LAP.get_trivial_solution([0.0, 0.0], my_t0, [sc_1_a_c * sc_1_a, 0.0], my_t1, init_guess = 2.0))

print("---- Trivial solution obtained ----")
my_LAP.trajectory_interpolate()


neighbors = []
N_neighbors = 10
max_deviation = 0.8

print("Initializing", N_neighbors, "neighbors with max. deviation =", max_deviation)
for i in range(N_neighbors):
    neighbors.append(LAP.copy(my_LAP))
    neighbors[-1].perturb_conditions([0.0, 0.0], [sc_1_a_c * sc_1_a , max_deviation / N_neighbors * (i + 1)])
    neighbors[-1].trajectory_interpolate()


print("---- Start the Metropolis - Hastings algorithm ----")

action_gradient_evolutions = []

for i in range(N_neighbors):
    print("Optimizing neighbor", i+1, "of", N_neighbors)
    action_gradient_evolutions.append(neighbors[i].metropolis_hastings(300, std = 0.01, bias = 3.0))

print("Printing the action gradient evolutions")
f_output = open("sc_1_a_g_evol.txt", "w")
for i_neigh in range(N_neighbors):
    for i_a_g in range(len(action_gradient_evolutions[i_neigh])):
        f_output.write("%f " % action_gradient_evolutions[i_neigh][i_a_g])
    f_output.write("\n")

f_output.close()

print("Obtaining spatial trajectory shapes")
f_spaces_output = open("sc_1_spaces.txt", "w")
my_spaces, my_timespace = my_LAP.get_path_geometry(100, True)

for i in range(len(my_timespace)):
    f_spaces_output.write("%f " % my_timespace[i])
f_spaces_output.write("\n")
for i in range(len(my_spaces[0])):
    f_spaces_output.write("%f," % my_spaces[0][i])
    f_spaces_output.write("%f " % my_spaces[1][i])
f_spaces_output.write("\n")

neighbors_spaces = []
for i in range(N_neighbors):
    neighbors_spaces.append(neighbors[i].get_path_geometry())
    for j in range(len(neighbors_spaces[-1][0])):
        f_spaces_output.write("%f," % neighbors_spaces[-1][0][j])
        f_spaces_output.write("%f " % neighbors_spaces[-1][1][j])
    f_spaces_output.write("\n")
f_spaces_output.close()



plt.title("Metropolis-Hastings, small monopole grid, 2 spatial dimensions\n (spline interpolation, 30 iterations)")

plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.plot(my_timespace, my_spaces[0], label = 'Trivial solution x')

for i in range(N_neighbors):
    cur_spaces = neighbors_spaces[i]
    plt.plot(my_timespace, cur_spaces[0], label = 'neigh. ' + str(i + 1) + ' x')
    plt.plot(my_timespace, cur_spaces[1], label = 'neigh. ' + str(i + 1) + ' y')

plt.legend(bbox_to_anchor=(1.05, 1))
plt.show()

plt.title("Spatial trajectories of the particles")
plt.plot(my_spaces[0], my_spaces[1], label = 'Trivial path')

for i in range(N_neighbors):
    cur_spaces = neighbors_spaces[i]
    plt.plot(cur_spaces[0], cur_spaces[1], label = 'neigh. ' + str(i + 1))

plt.legend()
plt.show()


end_dt = datetime.datetime.now()
duration = end_dt - start_dt
print("Process terminated in", duration.total_seconds(), "seconds.")


