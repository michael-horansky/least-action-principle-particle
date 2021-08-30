# least-action-principle-particle
Solving the problem of finding a particle's trajectory from the potential function using LAP and Monte Carlo

How to use this code:
  The file "action_integral_isotemporal.py" is the main piece of code. This is the solver of the LAP problem. It depeneds on the existence of the
  "monopole_potential_mesh.txt" file, which can be generated by running "monopole_potential_approximator.py". It generates two output files:
  "sc_1_a_g_evol.txt", which is a meta-analysis file showing how effective the process was, and "sc_1_spaces.txt", which are the actual trajectories
  of the investigated particles. These files should not be read, but rather plotted by running other short pieces of code.
  
  "spaces_plotter.py" can create plots from files like "sc_1_spaces.txt". It looks into a folder "evols/", and you have to rename the .txt file so that
  the number after "sc_" matches a variable inside the plotter code. Virtually the same goes for "a_g_evol_plotter.py", but this one creates plots from
  files like "sc_1_a_g_evol.txt". Again, place the text file into "evols/" folder and change the filename accordingly.
  
  Other code snippets: "potential_plotter.py" plots the potential mesh grid as a heatmap to visualize how the potential function looks.
  "trajectory_simulation.py" reads a "spaces" file and compares it with a simulation of a particle with the same starting dynamic properties
  using Newton's laws.

General notes:
  Everything should be in one folder, except for the text files which should be in the "evols/" folder within.
  The code requires python 3 and several libraries, which are listed below:
    numpy
    scipy
    matplotlib
    pynverse
    datetime

The code is provided as is, without license. It can be built upon, altered, or sold for profit. I don't really care.