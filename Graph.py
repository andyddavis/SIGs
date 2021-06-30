import numpy as np

from Domain import *
from Node import *

class Graph:                                                # Graph class

    def __init__(self, n, p_0, domain):
        self.n = n                                          # number of nodes on onee side of the graph
        self.p_0 = p_0                                      # probability of staying parameter
        self.domain = domain                                # domain (space coordinates and floe-velocity field)

        # initialize a list of nodes
        self.nodes = [Node(self.domain, k, self.n) for k in range(self.n**2)]
