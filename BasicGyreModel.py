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

    z = (x**2 + y**2)**0.5  # z is magnitude to normalise to unit vectors
    if (z == 0):                # catches division by zero
        z = 1
    return (y/z,-x/z)           # actual function (unit gyre, clockwise)


# O----------------------- define the graph  -----------------------O
class Node:                 # Node class
    # map from (i,j) in n x n to k in n^2 : (i,j) -> k= i*n + j
    # reverse map for k: k -> ((k-j)/n, k mod n)

    def __init__(self, k, n):
        if (k >= n**2):
            print("Warning: node index out of bounds.")
        self.n = n
        self.k = k
        self.j = k % n
        self.i = round((k-self.j)/n)

    # returns string representation including node labal k and indicies (i,j)
    def __str__(self):
        return("Node " + str(self.k) + "; Indices " + str((self.i,self.j)))

    # returns position of the node as a tuple (x,y); n-ary two-digit expression
    def pos(self):
        if (self.i >= self.n) or (self.j >= self.n):
            return print("Error 'pos_node': index out of bounds.")
        return( (self.j + 0.5) / self.n, 1 - (0.5 + self.i) / self.n)

    # returns the velocity vector at the position of the node
    def velocity(self):
        (x,y) = self.pos()
        (u,v) = floe_movement(x,y)
        return(u,v)

class Graph:                # Graph class
    # each node has a probability of mass staying (p_stay)
    # and a list of nodes equipped with a transition matrix

    def __init__(self, n, p_0 , mass):
        self.n = n                                          # number of nodes on onee side of the graph
        self.nodes = []                                     # list of nodes (not sure if needed)
        self.mass = mass_0                                  # initial coditions (mass at each node)
        self.p_0 = p_0                                      # probability of staying parameter
        self.tMatrix = np.zeros((n**2,n**2))                # transition matrix

    def initialise_nodes(self):
        self.nodes = []
        for k in range(0,self.n**2):
            self.nodes.append(Node(k,self.n))               # nodes are denoted with single index k

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
            if (k > (n-1)**2 and k < n**2 and v < 0):
                v = 0

            # pre-calculate useful values
            v_sum = abs(u) + abs(v)                         # sum of vector components
            v_mag = (u**2 + v**2)**0.5                      # magnitude of velocity vector at node
            p_stay = 1 / (1 + p_0 * v_mag)                  # probability of staying formula (credit: andy)
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

# O-------------------- simulation --------------------O
n = 10                              # number of nodes on one side (total nodes: n^2)
mass_0 = np.zeros(n**2)             # initial condition (inital mass at each node)
p_0 = 1                             # "inertial" parameter (changes probability of staying)

g = Graph(n, p_0, mass_0)           # create a graph
g.initialise_tMatrix()              # initialise it
