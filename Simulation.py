import numpy as np
import math

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

    def __init__(self, graph, mass, time, dt):
        self.graph = graph                                  # the graph
        self.mass = mass                                    # initial conditions (mass at each node)
        self.time = time                                    # total time of simulation
        self.dt = dt                                        # length of each timestep
        self.wind_tMatrix = np.zeros((graph.n**2,graph.n**2))           # wind transition matrix (external forcing)
        self.idle_tMatrix = np.zeros((graph.n**2,graph.n**2))           # diffusion transition matrix (internal forcing)

    # initialises the wind transition matrix (external forcing)
    def initialise_wind_tMatrix(self):
        boundary_condition = 0                              # (0) periodic or (1) closed
        for node in self.graph.nodes:                       # assign probabilities for eaech node:
            i = node.i; j = node.j; k = node.k              # gather node indices
            n = self.graph.n
            (u,v) = node.velocity()                         # gather component velocities
            v_mag = (u**2 + v**2)**0.5                      # magnitude of velocity vector at node

            # boundary conditions: closed (magnitude perserving - rotates velocity parallel to boundary)
            if boundary_condition == 1:
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
            p_i = self.graph.p_0 * n                             # p_stay parameter
            p_stay = 1 / (1 + p_i * v_mag)                  # probability of staying formula (credit: andy)
            self.wind_tMatrix[k,k] = p_stay                 # assign to diagonals

            # assign horizontal probabilities (adgacent to node k +/- 1 index)
            if (u > 0):
                if (j + 1 < n):
                    self.wind_tMatrix[k,k+1]   =  (1 - p_stay) * u / v_sum
                elif (boundary_condition == 0):
                    self.wind_tMatrix[k,k+1-n] =  (1 - p_stay) * u / v_sum
            elif (u < 0):
                if (j - 1 >= 0):
                    self.wind_tMatrix[k,k-1]   = -(1 - p_stay) * u / v_sum
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

    # initialise the diffusion matrix
    def initialise_idle_tMatrix(self):
        ## TODO: diffusion matrix
        pass

    # initialise the transition matrices
    def initialise_matrices(self):
        self.graph.initialise_nodes()           # initialise the nodes
        self.initialise_wind_tMatrix()          # initialise wind transition matrix
        self.initialise_idle_tMatrix()          # initialsie diffusion transition matrix

    def basic_sim(self):
        counter = 0
        frameskip = 5
        for timestep in range(0, round(self.time/self.dt)):
            # color scheme
            #cmap = 'gray'
            c_map = 'jet'

            # domain
            d = self.graph.domain
            dom = [d.x_lim[0], d.x_lim[1], d.y_lim[0], d.y_lim[1]]

            # plot
            if(counter % frameskip == 0):
                fig = plt.figure()
                plt.imshow(self.mass.reshape(self.graph.n,self.graph.n), cmap=c_map, interpolation='nearest',  extent=dom, vmin=0, vmax = 2 * np.sum(mass) / self.graph.n ** 2)
                plt.title("Time: " + str( timestep * self.dt))
                plt.colorbar()
                plt.savefig('figures/Step-'+str(int(timestep/frameskip)).zfill(10)+'.png', format='png', bbox_inches='tight')
                plt.close(fig)
            self.mass = self.mass.dot(self.graph.wind_tMatrix)
