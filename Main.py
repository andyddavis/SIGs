import numpy as np
import statistics as st
import math
import time as t
import matplotlib.pyplot as plt

from Domain import *
from Node import *
from Graph import *
from Simulation import *

from GyreDomain import *

# simulation attributes
options = dict()
options['timestep length'] = 1
options['total time'] = 100

# graph attributes
n = 50                               # number of nodes on one side (total nodes: n^2) 200
p_0 = 1  # "inertial" parameter (changes probability of staying)
domain = GyreDomain()
g = Graph(n, p_0, domain)       # create a graph

# check Courant number is maintained
C = 1                           # Courant number (?)
for node in g.nodes:
    for i in range(0,2):
        flag = False
        if (node.velocity()[i] * options['timestep length'] > C / n):
            print("Warning: dt not small enough: simulation may not reflect reality.")
            flag = True
            break
    if flag == True:
        break

# create and run simulation
sim = Simulation(g, options)    # create a GyreSimulation

sim.plot_steady_state()
#mass_0 = np.ones(n**2)
#sim.advection_sim(mass_0)
