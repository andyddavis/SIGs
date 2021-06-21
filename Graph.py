import numpy as np

from Domain import *
from Node import *

class Graph:                # Graph class
    # each node has a probability of mass staying (p_stay)
    # and a list of nodes equipped with a transition matrix

    def __init__(self, n, p_0 , mass):
        self.n = n                                          # number of nodes on onee side of the graph
        self.nodes = []                                     # list of nodes (not sure if needed)
        self.mass = mass                                    # initial coditions (mass at each node)
        self.p_0 = p_0                                      # probability of staying parameter
        self.tMatrix = np.zeros((n**2,n**2))                # transition matrix
        self.domain = Domain()

    def initialise_nodes(self):
        self.nodes = [None]*(self.n**2)
        for k in range(0,self.n**2):
            self.nodes[k] = Node(self.domain,k,self.n)               # nodes are denoted with single index k

    def initialise_tMatrix(self):
        self.initialise_nodes()                             # initialise the nodes

        for node in self.nodes:                             # assign probabilities for eaech node:

            k = node.k
            n = self.n
            (u,v) = node.velocity()                         # gather component velocities

            # boundary conditions: get rid of component if it points outside the space
            if (k % n == 0 and u < 0):
                u = 0
            if (k % n == (n-1) and u > 0):
                u = 0
            if (k < n and v > 0):
                v = 0
            if (k > (n-1)*n and k < n**2 and v < 0):
                v = 0

            # pre-calculate useful values
            v_sum = abs(u) + abs(v)                         # sum of vector components
            v_mag = (u**2 + v**2)**0.5                      # magnitude of velocity vector at node
            p_stay = 1 / (1 + self.p_0 * v_mag)                  # probability of staying formula (credit: andy)
            self.tMatrix[k,k] = p_stay            # assign to diagonals

            # assign horizontal probabilities (adgacent to node k +/- 1 index)
            if (u > 0):
                if(k % n + 1 < n):
                    self.tMatrix[k,k+1] = (1 - p_stay) * u / v_sum
            elif (u < 0):
                if(k % n - 1 >= 0):
                    self.tMatrix[k,k-1] = -(1 - p_stay) * u / v_sum
            # assign vertical probabilities (adgacent to node k +/- n index)
            if (v > 0):
                if(k - n >= 0):
                    self.tMatrix[k,k-n] = (1 - p_stay) * v / v_sum
            elif (v < 0):
                if(k + n < n ** 2):
                    self.tMatrix[k,k+n] = -(1 - p_stay) * v / v_sum
