import numpy as np
import statistics as st
import math

# O----------------- define the space (unit square)  -----------------O
x_lim = (0,1)
y_lim = (0,1)

# O---------- define floe movement function (vector field)  ----------O
def floe_movement(x,y):
    x = x-st.mean(x_lim)        # centers the graph at
    y = y-st.mean(y_lim)        # the center of limits

    z = math.sqrt(x**2 + y**2)  # z is magnitude to normalise to unit vectors
    if (z == 0):                # catches division by zero
        z = 1
    return (y/z,-x/z)           # actual function (unit gyre, clockwise)


# O----------------------- define the graph  -----------------------O
n = 3                       # number of nodes on one side of the
                            # ... square (total number of nodes: n^2)
p_stay = 1/4                # probability of staying at a node

class Node:                 # Node class
    # map from (i,j) in n x n to k in n^2 : (i,j) -> k= i*n + j
    # reverse map for k: k -> ((k-j)/n, k mod n)

    def __init__(self, k):
        if (k >= n**2):
            print("Warning: node index out of bounds.")
        self.k = k
        self.j = k % n
        self.i = round((k-self.j)/n)

    # returns string representation including node labal k and indicies (i,j)
    def __str__(self):
        return("Node " + str(self.k) + "; Indices " + str((self.i,self.j)))

    # returns position of the node as a tuple (x,y); n-ary two-digit expression
    def pos(self):
        if (self.i >= n) or (self.j >= n):
            return print("Error 'pos_node': index out of bounds.")
        return( (self.j + 0.5) / n, 1 - (0.5 + self.i) / n)

class Graph:                # Graph class
    # each node has a probability of mass staying (p_stay)
    # and a list of nodes equipped with a transition matrix

    def __init__(self, n, mass):
        self.n = n                                          # number of nodes on onee side of the graph
        self.p_stay = p_stay                                # probability of staying (not sure if needed)
        self.nodes = []                                     # list of nodes (not sure if needed)
        self.tMatrix = np.zeros((n**2,n**2))                # transition matrix
        self.mass = np.zeros(n**2)                          # initial coditions (mass at each node)

    def initialise_nodes(self):
        self.nodes = []
        for k in range(0,n**2):
            self.nodes.append(Node(k))                      # nodes are denoted with single index k

    def initialise_tMatrix(self):
        self.initialise_nodes()
        for node in self.nodes:
            self.tMatrix[node.k,node.k] = self.p_stay       # !!! BIG CHANGES NEEDED HERE !!!

# O-------------------- simulation --------------------O
mass = np.zeros(n**2)               # initial condition

g = Graph(n, mass)                 # create a graph
#g.initialise_tMatrix()              # initialise it
print(g.tMatrix)                    # check it
