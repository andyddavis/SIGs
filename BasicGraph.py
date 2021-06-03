import numpy as np
import time
import random as rand

# ~~~~~~~~~~~~~~~ a basic graph implementation ~~~~~~~~~~~~~~~ #

class BasicGraph:

    def __init__(self,size):
        self.size = size                        # size n for n nodes
        self.tMatrix = np.zeros((size,size))    # transition matrix

    def addEdge(self, node1, node2, p):
        if (node1 >= self.size) or (node2 >= self.size) or (abs(p)>1):
            print("Error: node/probability out of bounds.")
        else:
            self.tMatrix[node1][node2] = p      # assign edge with probability p


# ~~~~~~~~~~~~~~~ initialisation and testing ~~~~~~~~~~~~~~~ #

newGraph = BasicGraph(4)    # initalise the graph of size 4

for i in range(0,4):        # setting "stay" probabilities at 0.25
    newGraph.addEdge(i,i,0.25)

newGraph.addEdge(0,1,0.375)  # set up some edges manually
newGraph.addEdge(0,3,0.375)  # this matrix basically moves left and right with equal probability (3/8)

newGraph.addEdge(1,0,0.375)
newGraph.addEdge(1,2,0.375)

newGraph.addEdge(2,1,0.375)
newGraph.addEdge(2,3,0.375)

newGraph.addEdge(3,0,0.375)
newGraph.addEdge(3,2,0.375)

print(newGraph.tMatrix)     # view the transition matrix

# ~~~~~~~~~~~~~~~ printing where an agent is ~~~~~~~~~~~~~~~ #
def printMoves(agent):
    step = ""
    for i in range(0, newGraph.size):
        if i != agent:
            step += "o "
        else:
            step += "A "
    print(step)

# ~~~~~~~~~~~~~~~ random walk simulation ~~~~~~~~~~~~~~~ #
agent = 1                   # initialise start position of our lil dude
for i in range(0,10):
    val = rand.random()     # random value from 0 to 1
    counter = 0
    num = newGraph.tMatrix[agent][counter]
    while(num < val):       #cycle down the node's probability row til the right one is reached
        counter += 1
        num += newGraph.tMatrix[agent][counter]
    agent = counter         # move agent to the new spot
    printMoves(agent)
    time.sleep(0.5)
