import numpy as np
import statistics as st
import math
import time as t
import matplotlib.pyplot as plt

from Domain import *
from Node import *
from Graph import *
from Simulation import *

# time attributes
time = 100                      # time-length of simulation
dt = 1                          # size of each time-step

# graph attributes
n = 20                          # number of nodes on one side (total nodes: n^2) 200
mass_0=np.ones(n**2)            # initial condition (inital mass at each node)
p_0 = 1 * dt                        # "inertial" parameter (changes probability of staying)

g = Graph(n, p_0, mass_0)       # create a graph
g.initialise_graph()            # initialise it
#print(g.wind_tMatrix)

# check Courant number is maintained
C = 1                           # Courant number (?)
for node in g.nodes:
    if ((node.velocity()[0] * dt > C / n) or (node.velocity()[1] * dt > C / n)):
        print("Warning: 'dt' too large for proper simulation. Consider changing.")
        break

# create and run simulation
sim = Simulation(g, time, dt)    # create a GyreSimulation
sim.basic_sim()
