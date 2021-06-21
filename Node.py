import statistics as st

import weakref

from Domain import *

# O----------------------- define the graph  -----------------------O
class Node:                 # Node class
    # map from (i,j) in n x n to k in n^2 : (i,j) -> k= i*n + j
    # reverse map for k: k -> ((k-j)/n, k mod n)

    def __init__(self, domain, k, n):
        if (k >= n**2):
            print("Warning: node index out of bounds.")
        self.n = n
        self.k = k
        self.j = k % n
        self.i = round((k-self.j)/n)

        self.domain = weakref.ref(domain)

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
        (u,v) = self.domain().floe_movement(x,y)
        return(u,v)
