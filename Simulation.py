import numpy as np
import random as rand

from Advection import *
from Diffusion import *

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

class Simulation:

    def __init__(self, graph, options):
        self.graph = graph              # the graph (contains nodes and graph size/shape)
        self.options = options          # options for the simulation (total time, dt, frameskip)

    #--------------------------------------------------------------------------------------------------#

    # continuous advection with mass (repeated transition matrix application to initial mass condition)
    def advection_sim(self, mass):
        # initialise the advection
        advection = Advection(self.graph, self.options['timestep length'])
        # initialise the transition matrix
        advection.initialise_transition_matrix()

        # animation parameters <3
        # color scheme - 'gray', 'jet', etc
        c_map = 'jet'
        # domain
        d = self.graph.domain; dom = [d.x_lim[0], d.x_lim[1], d.y_lim[0], d.y_lim[1]]

        # simulation
        for timestep in range(0, round(self.options['total time']/self.options['timestep length'])):
            # plot every frameskip
            if(timestep % self.options['frameskip'] == 0):
                fig = plt.figure()
                plt.imshow(mass.reshape(self.graph.n,self.graph.n), cmap=c_map, interpolation='nearest',  extent=dom, vmin=0, vmax = 2 * np.sum(mass) / self.graph.n ** 2)
                plt.title("Time: " + str( timestep * self.options['timestep length']))
                plt.colorbar()
                plt.savefig('figures/Step-'+str(int(timestep/self.options['frameskip'])).zfill(10)+'.png', format='png', bbox_inches='tight')
                plt.close(fig)
            # evolve the mass each timestep
            mass = mass.dot(advection.transition_matrix)

    #--------------------------------------------------------------------------------------------------#

    # discrete advection with particles
    def advection_sim_particles(self):
        # initialise the advection
        advection = Advection(self.graph, self.options['timestep length'])
        # initialise the transition matrix
        advection.initialise_transition_matrix()

        # animation parameters <3
        # color scheme - 'gray', 'jet', etc
        c_map = 'gray'
        # domain
        d = self.graph.domain; dom = [d.x_lim[0], d.x_lim[1], d.y_lim[0], d.y_lim[1]]

        mass = np.zeros(self.graph.n**2) # mass array

        # simulation
        for timestep in range(0, round(self.options['total time']/self.options['timestep length'])):

            # calculate the "mass" at each node
            for i in range(0,self.graph.n**2):
                # just take the length of the particle list ( minus one for buffer None, assumes uniform particle mass of 1)
                mass[i] = len(self.graph.nodes[i].particles) - 1

            # plot every frameskip
            if(timestep % self.options['frameskip'] == 0):
                fig = plt.figure()
                plt.imshow(mass.reshape(self.graph.n,self.graph.n), cmap=c_map, interpolation='nearest',  extent=dom, vmin=0, vmax = np.sum(mass) / self.graph.n)
                plt.title("Time: " + str(timestep * self.options['timestep length']))
                plt.colorbar()
                plt.savefig('figures/Step-'+str(int(timestep/self.options['frameskip'])).zfill(10)+'.png', format='png', bbox_inches='tight')
                plt.close(fig)

            # simulating particle movement:
            # for each node, go through the particles and move them randomly according to transition matrix row probabilities
            for node in self.graph.nodes:
                # the None entry acts as a buffer to determine pre and post movement particles per timestep
                while node.particles[0] != None:
                    # choose random number [0,1] and move along the row (cnt) until you achieve the correct probability index
                    p = rand.random()
                    cnt = 0
                    dp = advection.transition_matrix[node.k, cnt]
                    while dp < p:
                        cnt += 1
                        dp += advection.transition_matrix[node.k, cnt]
                    # move the particle from the current node to the other (cnt) node, then remove it from the current node
                    self.graph.nodes[cnt].particles.append(node.particles[0])
                    del node.particles[0]
                # reset the buffer particle
                node.particles.remove(None); node.particles.append(None)

    #--------------------------------------------------------------------------------------------------#

    # continuous diffusion with mass (repeated transition matrix application to initial mass condition)
    def diffusion_sim(self, mass):
        # initialise the diffusion
        diffusion = Diffusion(self.graph, self.options['timestep length'], mass)

        # animation parameters <3
        # color scheme - 'gray', 'jet', etc
        c_map = 'jet'
        # domain
        d = self.graph.domain; dom = [d.x_lim[0], d.x_lim[1], d.y_lim[0], d.y_lim[1]]

        # simulation
        for timestep in range(0, round(self.options['total time']/self.options['timestep length'])):

            # plot every frameskip
            if(timestep % self.options['frameskip'] == 0):
                fig = plt.figure()
                plt.imshow(mass.reshape(self.graph.n,self.graph.n), cmap=c_map, interpolation='nearest',  extent=dom, vmin=0, vmax = 2 * np.sum(mass) / self.graph.n ** 2)
                plt.title("Time: " + str( timestep * self.options['timestep length']))
                plt.colorbar()
                plt.savefig('figures/Step-'+str(int(timestep/self.options['frameskip'])).zfill(10)+'.png', format='png', bbox_inches='tight')
                plt.close(fig)

            # update the diffusion's mass
            diffusion.mass = mass
            # update the diffusion transition matrix
            diffusion.initialise_transition_matrix()
            # evolve the mass each timestep
            mass = mass.dot(diffusion.transition_matrix)

    #--------------------------------------------------------------------------------------------------#
