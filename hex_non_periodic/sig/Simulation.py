import numpy as np
from scipy import interpolate

# import plotting packages and set default figure options
useserif = True # use a serif font with figures?
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
if useserif:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['text.usetex'] = False
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
        self.graph = Graph.create_hex_lattice(domain, options.get("num_nodes", 15))
        self.final_time = options.get("final_time", 1.0)
        self.num_timesteps = options.get("num_timesteps", 10)
        self.direc = options["direc"]
        self.nearest = [None]
        self.mass_evolution = []
        self.PE_evolution = []
        self.TE_evolution = []
        
        # # randomly choose where the particles will be placed
        # num_particles = options.get("num_particles", 1000)
        # placement = np.random.uniform(size=num_particles)
        # placement.sort()
        # node = 0
        # cumprob = self.graph.nodes[0].mass_density
        # for place in placement:
        #     if place>cumprob:
        #         # recompute the mass density using the particle estimte
        #         self.graph.nodes[node].mass_density = self.graph.nodes[node].compute_mass_density(num_particles)
        #         self.graph.nodes[node].update_energy()

        #         node += 1
        #         cumprob += self.graph.nodes[node].mass_density

        #     # create a particle to add to the current node
        #     self.graph.nodes[node].add_particle(np.random.normal(size=2), 0)

        # # recompute the mass density using the particle estimate
        # for i in range(node-1, len(self.graph.nodes)):
        #     self.graph.nodes[i].mass_density = self.graph.nodes[i].compute_mass_density(num_particles)
        #     self.graph.nodes[i].update_energy()
        
        # 8 neighbors center diffusion:
        
        
        # num_particles = options.get("num_particles", 1000)
        # num_nodes = options.get("num_nodes", 15)
        # for i in range(num_nodes*num_nodes):
        #     if( i == 4949 or i == 4950 or i==5049 or i == 5050 ):    # if( i == 779 or i == 780 or i==819 or i == 820 ):   # i == 4949 or i == 4950 or i==5049 or i == 5050 
        #         for j in range(int(num_particles/4)):
        #             self.graph.nodes[i].add_particle(np.random.normal(size=2), 0)

        num_particles = options.get("num_particles", 1000)
        num_nodes = options.get("num_nodes", 15)
        
        
#    ++++++++++++++++++++++++++++++++++++++
#    initialize all particles at the center of the graph
        for i in range(num_particles):
            self.graph.nodes[num_nodes][num_nodes].add_particle(np.random.normal(size=2),0)   #  np.random.normal(size=2), 0
#    ++++++++++++++++++++++++++++++++++++++

#    ++++++++++++++++++++++++++++++++++++++
#    same number of particles at every node 
        # for node in self.graph.nodes.flatten():
        #     if(node):
        #         for i in range(10):
        #             node.add_particle(np.random.normal(size=2), 0)
#    ++++++++++++++++++++++++++++++++++++++

#  compute energy and mass density
        for node in self.graph.nodes.flatten():
            if (node):
                node.mass_density = node.compute_mass_density(num_particles)
                node.update_energy()
        # for i in range(num_nodes*num_nodes):
        #     self.graph.nodes[i].mass_density = self.graph.nodes[i].compute_mass_density(num_particles)
        #     self.graph.nodes[i].update_energy()
        self.initial_PE = self.get_initial_PE()
    def run(self):
        dt = self.final_time/self.num_timesteps
        time = 0.0
        zfill = int(np.log10(self.num_timesteps)+2)
        #self.domain.plot_external_velocity(0, dt, zfill)
        self.plot_mass_density(filename='mass_density-'+str(0).zfill(zfill), ext='png')
        self.plot_PE(filename='potential_energy-'+str(0).zfill(zfill), ext='png')

        for t in range(self.num_timesteps):
            print(str(t).zfill(zfill))
            time += dt
            #self.graph.acceleration_collisions(dt)
            self.graph.new_move_particles(dt, self.initial_PE)#move_particles(dt, self.initial_PE)
            #self.graph.increase_time(dt)
            #self.mass_evolution.append(self.graph.get_mass())
            #self.PE_evolution.append(self.graph.get_PE())
            #self.TE_evolution.append(self.graph.get_TE())
            #self.domain.plot_external_velocity(self.graph.nodes[num_nodes][num_nodes].time, dt, zfill)
            self.plot_mass_density(filename='mass_density-'+str(t+1).zfill(zfill), ext='png')
            # self.plot_potential_energy(filename='potential_energy-'+str(t+1).zfill(zfill), ext='png')
        print(time)


    def plot_mass_density(self, filename='mass_density', ext='pdf', resolution=100):
        center_x, center_y = self.graph.compute_center_position()
        PE_x, PE_y = self.graph.compute_center_PE()
        data=self.graph.create_scatter_data()
        fig=plt.figure(figsize=(10,10))
        plt.scatter(data[0], data[1],alpha=0.3,cmap='coolwarm')
        plt.scatter(center_x,center_y, s=50, c='red',
                    label=f'mass center: x= {center_x}, y={center_y}')
        plt.scatter(PE_x,PE_y, s=50, c='green',
                    label=f'\nPE center: x= {PE_x}, y={PE_y}')
        plt.scatter([],[],)
        plt.xlim(-1.0, 1.0)
        plt.ylim(-1.0, 1.0)
        plt.legend()
        fig.savefig(self.direc + filename, dpi=fig.dpi)
        plt.close()
        
        
    def plot_PE(self, filename='PE', ext='pdf', resolution=100):
        center_x, center_y = self.graph.compute_center_PE()
        PE_x, PE_y = self.graph.compute_center_PE()
        data=self.graph.create_scatter_data()
        fig=plt.figure(figsize=(10,10))
        plt.scatter(data[0], data[1],alpha=0.3,cmap='coolwarm')
        plt.scatter(center_x,center_y, s=50, c='red',
                    label=f'mass center: x= {center_x}, y={center_y}')
        plt.scatter(PE_x,PE_y, s=50, c='green',
                    label=f'\nPE center: x= {PE_x}, y={PE_y}')
        plt.scatter([],[],)
        plt.xlim(-1.0, 1.0)
        plt.ylim(-1.0, 1.0)
        plt.legend()
        fig.savefig(self.direc + filename, dpi=fig.dpi)
        plt.close()
        
        
    

    def plot_potential_energy(self, filename='potential_energy', ext='pdf', resolution=100):
        zdata = [0.0]*(len(self.graph.nodes)**2)
        points = [None]*(len(self.graph.nodes)**2)
        for i in range(len(self.graph.nodes.flatten())):
            node = self.graph.nodes.flatten()[i]
            if(node):
                points[i] = node.position
                zdata[i] = node.potential_energy
            else:
                continue
        X, Y = np.meshgrid(np.linspace(-self.domain.radius, self.domain.radius, resolution), 
                           np.linspace(-self.domain.radius, self.domain.radius, resolution))
        center_x, center_y = self.graph.compute_center_position()
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
        plt.savefig(self.direc + filename+'.'+ext, format=ext, bbox_inches='tight')
        plt.close(fig)

    def get_initial_PE (self):
        counter = 0
        for node in self.graph.nodes.flatten():
            if (node):
                counter += len(node.particles)*node.potential_energy
        return counter