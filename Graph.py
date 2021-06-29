import numpy as np

from Domain import *
from Node import *

class Graph:                                                # Graph class

    def __init__(self, n, p_0, domain):
        self.n = n                                          # number of nodes on onee side of the graph
        self.nodes = []                                     # list of nodes (not sure if needed)
        #self.mass = mass                                    # initial conditions (mass at each node)
        self.p_0 = p_0                                      # probability of staying parameter
        #self.wind_tMatrix = np.zeros((n**2,n**2))           # wind transition matrix (external forcing)
        #self.idle_tMatrix = np.zeros((n**2,n**2))           # diffusion transition matrix (internal forcing)
        self.domain = domain                                # domain (space coordinates and floe-velocity field)

    # creates the list of nodes in the graph
    def initialise_nodes(self):
        self.nodes = [None]*(self.n**2)
        for k in range(0,self.n**2):
            self.nodes[k] = Node(self.domain,k,self.n)      # nodes are denoted with single index k

    # initialises the wind transition matrix
    def initialise_wind_tMatrix(self):
        boundary_condition = 0                              # (0) periodic or (1) closed
        for node in self.nodes:                             # assign probabilities for eaech node:
            i = node.i; j = node.j; k = node.k              # gather node indices
            n = self.n
            (u,v) = node.velocity()                         # gather component velocities
            v_mag = (u**2 + v**2)**0.5                      # magnitude of velocity vector at node

            # boundary conditions: closed (magnitude perserving)
            if boundary_condition == 0:
                if (j == 0 and u < 0):
                    u = 0; v = np.sign(v) * v_mag
                if (j == (n-1) and u > 0):
                    u = 0; v = np.sign(v) * v_mag
                if (i == 0 and v > 0):
                    v = 0; u = np.sign(u) * v_mag
                if (i == (n-1) and v < 0):
                    v = 0; u = np.sign(u) * v_mag

            # pre-calculate useful values
            v_sum = abs(u) + abs(v)                         # sum of vector components
            p_i = self.p_0 * n                              # p_stay parameter
            p_stay = 1 / (1 + p_i * v_mag)                  # probability of staying formula (credit: andy)
            self.wind_tMatrix[k,k] = p_stay                 # assign to diagonals

            # assign horizontal probabilities (adgacent to node k +/- 1 index)
            if (u > 0):
                if (j + 1 < n):
                    self.wind_tMatrix[k,k+1] = (1 - p_stay) * u / v_sum
                elif (boundary_condition == 0):
                    self.wind_tMatrix[k,k+1-n] = (1 - p_stay) * u / v_sum
            elif (u < 0):
                if (j - 1 >= 0):
                    self.wind_tMatrix[k,k-1] = -(1 - p_stay) * u / v_sum
                elif (boundary_condition == 0):
                    self.wind_tMatrix[k,k-1+n] = -(1 - p_stay) * u / v_sum
            # assign vertical probabilities (adgacent to node k +/- n index)
            if (v > 0):
                if (i - 1 >= 0):
                    self.wind_tMatrix[k,k-n] = (1 - p_stay) * v / v_sum
                elif (boundary_condition == 0):
                    self.wind_tMatrix[k,k-n+n**2] = (1 - p_stay) * v / v_sum
            elif (v < 0):
                if (i + 1 < n):
                    self.wind_tMatrix[k,k+n] = -(1 - p_stay) * v / v_sum
                elif (boundary_condition == 0):
                    self.wind_tMatrix[k,k+n-n**2] = -(1 - p_stay) * v / v_sum

    def initialise_idle_tMatrix(self):
        ## TODO: diffusion matrix
        pass

    def initialise_graph(self):
        self.initialise_nodes()                 # initialise the nodes
        self.initialise_wind_tMatrix()          # initialise wind transition matrix
        self.initialise_idle_tMatrix()          # initialsie diffusion transition matrix
