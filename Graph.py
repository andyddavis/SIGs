import numpy as np

from Domain import *
from Node import *

class Graph:                                                # Graph class

    def __init__(self, n, p_0, domain):
        self.n = n                                          # number of nodes on onee side of the graph
        self.nodes = []                                     # list of nodes
        self.p_0 = p_0                                      # probability of staying parameter
        self.domain = domain                                # domain (space coordinates and floe-velocity field)

    # creates the list of nodes in the graph
    def initialise_nodes(self):
        self.nodes = [None]*(self.n**2)
        for k in range(0,self.n**2):
            self.nodes[k] = Node(self.domain, k, self.n)      # nodes are denoted with single index k
