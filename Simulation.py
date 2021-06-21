import numpy as np
import statistics as st
import math
import time as t

# import plotting packages and set default figure options
useserif = False # use a serif font with figures?
import matplotlib as mpl
import matplotlib.pyplot as plt
if useserif:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['text.usetex'] = True
else:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams['text.usetex'] = False
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

# O-------------------- simulation --------------------O
class Simulation:

    def __init__(self, graph, time, dt):
        self.graph = graph
        self.time = time
        self.dt = dt

    def basic_sim(self):
        m = self.graph.mass
        for timestep in range(0, round((self.time+self.dt)/self.dt)):
            #print("time: ", timestep * dt)
            #print(m.reshape(n,n))
            #plt.imshow(m.reshape(n,n), cmap='gray', interpolation='nearest',  extent=[0, 1, 0, 1], vmin=0, vmax=1)
            #plt.colorbar()
            #plt.title("Time: " + str( timestep * dt))
            #plt.draw()
            #plt.pause(self.dt)
            m = m.dot(self.graph.tMatrix)
            #t.sleep(self.dt)
            # plt.show()

            fig = plt.figure()
            plt.imshow(m.reshape(self.graph.n,self.graph.n), cmap='jet', interpolation='nearest',  extent=[0, 1, 0, 1])
            plt.colorbar()
            plt.savefig('figures/Step-'+str(timestep).zfill(10)+'.png', format='png', bbox_inches='tight')
            plt.close(fig)
