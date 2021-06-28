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
            # color scheme
            #cmap = 'gray'
            c_map = 'jet'

            # domain
            d = self.graph.domain
            dom = [d.x_lim[0], d.x_lim[1], d.y_lim[0], d.y_lim[1]]

            # plot
            fig = plt.figure()
            plt.imshow(m.reshape(self.graph.n,self.graph.n), cmap=c_map, interpolation='nearest',  extent=dom, vmin=0, vmax = 2 * np.sum(m) / self.graph.n ** 2)
            plt.title("Time: " + str( timestep * self.dt))
            plt.colorbar()
            plt.savefig('figures/Step-'+str(timestep).zfill(10)+'.png', format='png', bbox_inches='tight')
            plt.close(fig)
            m = m.dot(self.graph.wind_tMatrix)
