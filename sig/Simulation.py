import numpy as np
from scipy import interpolate

# import plotting packages and set default figure options
useserif = True # use a serif font with figures?
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
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

from sig.Graph import *

class Simulation:
    def __init__(self, domain, options):
        self.domain = domain
        self.graph = Graph.create_square_lattice(domain, options.get("num_nodes", 15))
        self.final_time = options.get("final_time", 1.0)
        self.num_timesteps = options.get("num_timesteps", 10)
        self.nearest = [None]

        # randomly choose where the particles will be placed
        num_particles = options.get("num_particles", 1000)
        placement = np.random.uniform(size=num_particles)
        placement.sort()
        node = 0
        cumprob = self.graph.nodes[0].mass_density
        for place in placement:
            if place>cumprob:
                # recompute the mass density using the particle estimte
                self.graph.nodes[node].mass_density = self.graph.nodes[node].compute_mass_density(num_particles)
                self.graph.nodes[node].update_energy()

                node += 1
                cumprob += self.graph.nodes[node].mass_density

            # create a particle to add to the current node
            self.graph.nodes[node].add_particle(np.random.normal(size=2))

        # recompute the mass density using the particle estimate
        for i in range(node-1, len(self.graph.nodes)):
            self.graph.nodes[i].mass_density = self.graph.nodes[i].compute_mass_density(num_particles)
            self.graph.nodes[i].update_energy()

    def run(self):
        dt = self.final_time/self.num_timesteps
        time = 0.0
        zfill = int(np.log10(self.num_timesteps)+2)

        self.plot_mass_density(filename='periodic_gyres_mass/mass_density-'+str(0).zfill(zfill), ext='png')
        self.plot_potential_energy(filename='periodic_gyres_pe/potential_energy-'+str(0).zfill(zfill), ext='png')

        for t in range(self.num_timesteps):
            print(str(t).zfill(zfill))
            time += dt
            self.graph.acceleration_collision(dt)

            self.graph.move_particles(dt)

            self.plot_mass_density(filename='periodic_gyres_mass/mass_density-'+str(t+1).zfill(zfill), ext='png')
            self.plot_potential_energy(filename='periodic_gyres_pe/potential_energy-'+str(t+1).zfill(zfill), ext='png')
        print(time)

    def plot_mass_density(self, filename='mass_density', ext='pdf', resolution=100):
        zdata = [0.0]*len(self.graph.nodes)
        points = [None]*len(self.graph.nodes)
        for i in range(len(self.graph.nodes)):
            points[i] = self.graph.nodes[i].position
            zdata[i] = self.graph.nodes[i].mass_density

        X, Y = np.meshgrid(np.linspace(self.domain.x_lim[0], self.domain.x_lim[1], resolution), np.linspace(self.domain.y_lim[0], self.domain.y_lim[1], resolution))
        mass_density = interpolate.griddata(points, zdata, (X, Y), method='nearest')

        fig = plt.figure()
        ax = plt.gca()
        pc = ax.pcolor(X, Y, mass_density, cmap='jet')
        fig.colorbar(pc)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel(r'$x$')
        ax.set_xlabel(r'$y$')
        plt.savefig(filename+'.'+ext, format=ext, bbox_inches='tight')
        plt.close(fig)

    def plot_potential_energy(self, filename='potential_energy', ext='pdf', resolution=100):
        zdata = [0.0]*len(self.graph.nodes)
        points = [None]*len(self.graph.nodes)
        for i in range(len(self.graph.nodes)):
            points[i] = self.graph.nodes[i].position
            zdata[i] = self.graph.nodes[i].potential_energy

        X, Y = np.meshgrid(np.linspace(self.domain.x_lim[0], self.domain.x_lim[1], resolution), np.linspace(self.domain.y_lim[0], self.domain.y_lim[1], resolution))
        potential_energy = interpolate.griddata(points, zdata, (X, Y), method='nearest')

        fig = plt.figure()
        ax = plt.gca()
        pc = ax.pcolor(X, Y, potential_energy, cmap='jet')
        fig.colorbar(pc)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel(r'$x$')
        ax.set_xlabel(r'$y$')
        plt.savefig(filename+'.'+ext, format=ext, bbox_inches='tight')
        plt.close(fig)
