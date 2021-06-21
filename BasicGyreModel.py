import numpy as np
import statistics as st
import math
import time as t
import matplotlib.pyplot as plt

from Domain import *
from Node import *
from Graph import *
from Simulation import *

time = 1000                           # time-length of simulation
dt = 1                               # size of each time-step

n = 20                               # number of nodes on one side (total nodes: n^2)
#mass_0=np.zeros(n**2); mass_0[round(n**2 / 3)]= 2*n  # initial condition (inital mass at each node)
mass_0=np.ones(n**2); # initial condition (inital mass at each node)
p_0 = dt * 2                         # "inertial" parameter (changes probability of staying)

g = Graph(n, p_0, mass_0)            # create a graph
g.initialise_tMatrix()               # initialise it
#print(g.tMatrix)

sim = Simulation(g, time, dt)    # create a GyreSimulation
sim.basic_sim()
