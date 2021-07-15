import numpy as np
import statistics as st
import math
import time as t
import matplotlib.pyplot as plt

from Domain import *
from GyreDomain import *
from Particle import *
from Node import *
from Graph import *
from Simulation import *


# simulation attributes
options = dict()
options['timestep length'] = 1
options['total time'] = 1000
options['frameskip'] = 5

# graph attributes
n = 50                  # number of nodes on one side (total nodes: n^2) 200
p_0 = 1                 # "inertial" parameter (changes probability of staying)
domain = GyreDomain()   # domain of choice
g = Graph(n, p_0, domain)

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

# particle placement for particle sim
for i in range(0,300):
    g.nodes[int(n**2 / 3 + round(2*n/3))].particles.insert(0,Particle(0,1))

# mass distribution for mass sim
#mass_0 = np.ones(n**2)
mass_0 = np.zeros(n**2)
mass_0[int((n**2 + n )/ 2)] = 100

# create and run simulation
sim = Simulation(g, options)    # create a GyreSimulation

#domain.plot()
#sim.advection_sim(mass_0)
#sim.diffusion_sim(mass_0)
sim.advection_sim_particles()
