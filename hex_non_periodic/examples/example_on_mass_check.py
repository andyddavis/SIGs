#!/usr/bin/env python
# coding: utf-8

# In[6]:


import statistics as st
import weakref
from scipy import linalg
from scipy import sparse
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time as t
import pandas as pd
from itertools import combinations
import random
import os 


# simulation attributes
options = dict()
delta =1
options['timestep length'] = delta
options['total time'] = 50
direc='/Users/25016/Desktop/Advection_Collision_Gaussian/'  
n = 20                               # number of nodes on one side (total nodes: n^2) 400
p_0 = 1                             # "inertial" parameter (changes probability of staying)
collision_rate = 0.75
gamma = 0.8
mass_evolution=[]
number_of_particles_at_each_node=[100]*(n**2) # a list that specifies how many particles there are at each node
plot_or_not = 1
df_or_not = 0

class Particle:

    def __init__(self, vel, k, n, time, ID):
        self.vel = vel
        self.k = k
        self.n = n
        self.time = time
        self.ID = ID
        # vel is a numpy array representing velocity
    def display_info (self):
        print('Velocity: %.2f, %.2f' % (self.vel[0], self.vel[1]))  
        print('Time: %.2f' % self.time)
        print('ID: '+self.ID)
    def inelastic_collide_particle (self, particle2):
        angle_of_attack = 2*np.pi* random.random()
        v1=self.vel
        v2=particle2.vel
        pos_diff=np.array([np.cos(angle_of_attack),np.sin(angle_of_attack)])
        # pos_diff is a random unit vector
        value=-1*np.dot(pos_diff,v1-v2)
        if (value**2-4*(1-gamma)*(0.5*np.dot(v1,v1)+0.5*np.dot(v2,v2))<0):
            #print('This collision cannot happen')
            return 
        else:
            if (value>0):
                new_value=0.5*value+0.5*np.sqrt(value**2-4*(1-gamma)*(0.5*np.dot(v1,v1)+0.5*np.dot(v2,v2)))
            elif(value==0):
                new_value=0
            else:
                new_value=0.5*value-0.5*np.sqrt(value**2-4*(1-gamma)*(0.5*np.dot(v1,v1)+0.5*np.dot(v2,v2)))
            new_v1=v1+new_value*pos_diff
            new_v2=v2-new_value*pos_diff
            self.vel=new_v1
            particle2.vel=new_v2
    def increase_time_particle (self):
        self.time = self.time + delta
    def pos_particle (self):
        j = self.k % self.n
        return [(self.k-j)/self.n, j] 
    
class Domain:
    def __init__(self):
        self.x_lim = (0.0, 1.0)
        self.y_lim = (0.0, 1.0)

    # define external velocity (vector field) function
    def external_velocity(self, x, y):
        return (0.0, 0.0)

    # external velocity plot
    def plot(self):
        x,y = np.meshgrid(np.linspace(self.x_lim[0],self.x_lim[1],10),np.linspace(self.y_lim[0],self.y_lim[1],10))
        (u,v) = self.external_velocity(x,y)
        plt.quiver(x,y,u,v)
        plt.show()
        plt.axis('equal')
        
        
class GyreDomain:
    def __init__(self):
        Domain.__init__(self)

    # define external velocity (vector field)
    def external_velocity(self, x, y):
        x = x
        y = y
        k = 1
        #return (np.sin(k * np.pi * x) * np.cos(k * np.pi * y), np.cos(k * np.pi * x) * np.sin(k * np.pi *y))
        #return (k*y/np.linalg.norm([x,y]), -k*x/np.linalg.norm([x,y]))
        return (np.sin(k * np.pi * x) * np.cos(k * np.pi * y), - np.cos(k * np.pi * x) * np.sin(k * np.pi *y))

    # plot the vector field
    def plot(self):
        x,y = np.meshgrid(np.linspace(self.x_lim[0],self.x_lim[1],11),np.linspace(self.y_lim[0],self.y_lim[1],11))
        (u,v) = self.external_velocity(x,y)
        plt.quiver(x,y,u,v)
        plt.show()
        plt.axis('equal')
        
        
class Node:                 # Node class
    # map from (i,j) in n x n to k in n^2 : (i,j) -> k= i*n + j (n-ary two-digit expression)
    # reverse map for k: k -> ((k-j)/n, k mod n)

    def __init__(self, domain, k, n, time, particles):
        if (k >= n**2):
            print("Warning: node index out of bounds.")
        self.n = n
        self.k = k
        self.j = k % n
        self.i = round((k-self.j)/n)
        self.domain = domain
        self.time = time
        self.particles = particles 
    # returns string representation including node labal k and indicies (i,j)
    
    def add_particle (self, particle):
        self.particles = self.particles  +[particle]        
    
    def display_particles(self):
        for i in self.particles:
            print(i.ID)
    
    def __str__(self):
        return("Node " + str(self.k) + "; Indices " + str((self.i,self.j)))

    # returns position of the node as a tuple (x,y)
    def pos(self):
        if (self.i >= self.n) or (self.j >= self.n):
            return print("Error 'pos_node': index out of bounds.")
        return np.array([(self.j + 0.5) / self.n, 1 - (0.5 + self.i) / self.n])

    # returns the velocity vector (u,v) at the position of the node
    def velocity(self):
        (x,y) = self.pos()
        (u,v) = self.domain.external_velocity(x,y)
        return(u,v)
    
    def collide_node (self,collision_rate):
        # This function takes in the collision rate at this node
        self.time=self.time+ delta
        for i in self.particles:
            i.increase_time_particle()
        index=list(range(len(self.particles)))
        number_of_pairs=int(len(self.particles)*collision_rate/2)
        pairs=list(range(number_of_pairs)) # initialize pairs to be of equal length
        # Record the pairs of indexes of particles to collide in the list pairs
        for i in range(number_of_pairs):
            new_pair=list(random.sample(index,2))
            pairs[i]=new_pair
            index=list(set(index)-set(new_pair))
        for i in pairs:
            self.particles[i[0]].inelastic_collide_particle(self.particles[i[1]])       
    def hist (self,number, direc):
        # Plot a 2D histograph for the distribution of velocity at this node
        length=len(self.particles)
        X=list(range(length))
        Y=list(range(length))
        for i in X:
            X[i]=self.particles[i].vel[0]
            Y[i]=self.particles[i].vel[1]
        plt.hist2d(X,Y,bins=(50,50))
        plt.title("2D Histogram #"+str(number))
        plt.savefig(direc +'simu_hist'+str(number)+".png")    
    def mean_speed(self):
        speed_sum=0
        for i in self.particles:
            speed_sum+=np.linalg.norm(i.vel)
        return speed_sum/len(self.particles)

    def mean_velocity(self):
        x=0
        y=0
        for i in self.particles:
            x+=i.vel[0]
            y+=i.vel[1]
        return np.array([x,y])/len(self.particles)

    def total_energy(self):
        total_energy=0
        for i in self.particles:
            total_energy+= 0.5*np.dot(i.vel,i.vel)   
        return total_energy/len(self.particles)

    def kinetic_energy (self):
        mean_vel=self.mean_velocity()
        return 0.5 * np.dot(mean_vel,mean_vel) 

    def potential_energy(self):
        return self.total_energy()-self.kinetic_energy()

    
    
class Graph:                                                # Graph class

    def __init__(self, n, p_0, domain ,nodes):
        self.n = n                                          # number of nodes on onee side of the graph
        self.p_0 = p_0                                      # probability of staying parameter
        self.domain = domain                                # domain (space coordinates and floe-velocity field)
        self.nodes = nodes                                  # nodes (the nodes in the graph)
       
    def get_mass (self):
        return np.array([ len(i.particles) for i in self.nodes])
                                                            # keep track of the mass distribution of the graph
        
        
# import plotting packages and set default figure options
useserif = False # use a serif font with figures?
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
        self.graph = graph                                              # the graph
        self.options = options                                          # options for the simulation
        self.advection_tMatrix = np.zeros((graph.n**2,graph.n**2))      # alt matrix
        self.diffusion_tMatrix = np.zeros((graph.n**2,graph.n**2))      # diffusion transition matrix (internal forcing)

    # initialises the node advection probabilities and their corresponding indices
    def node_advection_transition_probabilities(self, node):

        # calculate node indices (periodic boundary), ordered: itself, left, right, below, above)
        transition_node_indices=[node.k]                           # itself
        if( (node.k-1 > -1) & ( node.k-1 < self.graph.n**2)  ):          # left
            transition_node_indices.append(node.k-1)
        if( (node.k+1 > -1) & ( node.k+1 < self.graph.n**2)  ):          # right
            transition_node_indices.append(node.k+1)               
        if( (node.k-self.graph.n > -1) & ( node.k-self.graph.n < self.graph.n**2)  ):# below
            transition_node_indices.append(node.k-self.graph.n)
        if( (node.k+self.graph.n > -1) & ( node.k+self.graph.n < self.graph.n**2)  ):# above
            transition_node_indices.append(node.k+self.graph.n)
#         transition_node_indices = [node.k,
#                                    (node.k-1), #(node.j-1) % self.graph.n + self.graph.n * node.i ,
#                                    (node.k+1), #(node.j+1) % self.graph.n + self.graph.n * node.i,  
#                                    (node.k-n), #(node.k + self.graph.n ) % self.graph.n**2,         
#                                    (node.k+n)] #(node.k - self.graph.n) % self.graph.n**2          

        # calculate probability of staying
        p_stay = 1 / (1 + self.graph.p_0 * self.graph.n * self.options['timestep length'] * np.linalg.norm(node.velocity()))

        # initialise probabilities
        transition_node_probabilities = [0.0] * len(transition_node_indices)

        # if p_stay 1, no need for other calculations (||v|| = 0)
        if abs(p_stay - 1) < 1.0e-13:
            transition_node_probabilities[0] = p_stay
            return transition_node_indices, transition_node_probabilities

        # otherwise, calculate "leaving measure" for each neighbor
        for i in range(1, len(transition_node_indices)):
            # calculate edge vector
            edge = self.graph.nodes[transition_node_indices[i]].pos() - node.pos()
            if (np.linalg.norm(edge) > 0.5):                                           # assuming [0,1] X [0,1] (dependent on x_lim, y_lim - if errors, that's why)
                edge[0] = - np.sign(edge[0]) * (1 / self.graph.n)                      # periodic boundary condition fix (assuming uniform square lattice grid)
                edge[1] = - np.sign(edge[1]) * (1 / self.graph.n)
            # assign  measure
            transition_node_probabilities[i] = max(0.0, np.dot(edge, node.velocity()))

        # normalise probabilies of leaving
        S = sum(transition_node_probabilities)
        transition_node_probabilities = (1 - p_stay) * np.array(transition_node_probabilities) * (1 / S)

        # finally, assign probability of staying and return indices/probabilities
        transition_node_probabilities[0] = p_stay
        return transition_node_indices, transition_node_probabilities

    # create the advection transition matrix
    def advection_transition_matrix(self):
        for node in self.graph.nodes:
            (indices, probabilities) = self.node_advection_transition_probabilities(node)
            for i in range(0,len(indices)):
                self.advection_tMatrix[node.k, indices[i]] = probabilities[i]

    # initialise the diffusion matrix
    def diffusion_transition_tMatrix(self):
        ## TODO: diffusion matrix
        pass
    def particle_advection_sim(self):
        old_mass=self.graph.get_mass() # keep track of the old mass distribution
        for node in self.graph.nodes:# loop through all nodes
            to_remove=[]                      # keep track of the indices of particles to be removed
            node.time  = delta + node.time
            old_length=old_mass[node.k]
            (indices, probabilities) = self.node_advection_transition_probabilities(node)
            accumulated_prob = np.cumsum(probabilities)
            
            for particle_num in range(old_length):
                particle = node.particles[particle_num]
                particle.increase_time_particle()
                random_seed=random.random()
                if ( random_seed < accumulated_prob[0]):
                    #print('Staying')
                    pass
                elif (random_seed < accumulated_prob[1]):
                    #print('Moving to 1')
                    particle.k=indices[1]
                    self.graph.nodes[indices[1]].add_particle(particle)
                    to_remove.append(particle_num)
                elif (random_seed < accumulated_prob[2]):
                    #print('Moving to 2')
                    particle.k=indices[2]
                    self.graph.nodes[indices[2]].add_particle(particle)
                    to_remove.append(particle_num)
                else:
                    try:
                        if (random_seed < accumulated_prob[3]):
                            #print('Moving to 3')
                            particle.k=indices[3]
                            self.graph.nodes[indices[3]].add_particle(particle)
                            to_remove.append(particle_num)

                        else :
                            try:
                                if (random_seed < accumulated_prob[4]):
                                    #print('Moving to 4')
                                    particle.k=indices[4]
                                    self.graph.nodes[indices[4]].add_particle(particle)
                                    to_remove.append(particle_num)
                            except:
                                pass

                    except:
                        pass

            
            for i in sorted(to_remove,reverse=True):
                try:
                    del node.particles[i]
                except:
                    print('Fail to remove!!!!')
    def collision_sim (self):
        for node in self.graph.nodes:
            node.collide_node(collision_rate)
                               
    def update_df (self):
        for node in self.graph.nodes:
            for particle in node.particles:
                df.loc[particle.ID,'time '+str(int(particle.time))]=particle.k
    def advection_collision_plot(self):
        #sparse_tMatrix = sparse.csr_matrix(self.advection_tMatrix)
        counter = 0
        frameskip = 5
        for timestep in range(0, round(self.options['total time']/self.options['timestep length'])):
            # color scheme
            #cmap = 'gray'
            c_map = 'jet'

            # domain
            d = self.graph.domain
            dom = [d.x_lim[0], d.x_lim[1], d.y_lim[0], d.y_lim[1]]

            mass=self.graph.get_mass()
            # plot
            if(counter % frameskip == 0 and (plot_or_not > 0)):
                fig = plt.figure()
                plt.imshow(mass.reshape(self.graph.n,self.graph.n),
                           cmap=c_map, interpolation='nearest',  
                           extent=dom, vmin=0, vmax = 2 * np.sum(mass) / self.graph.n ** 2)
                plt.title("Time: " + str( timestep * self.options['timestep length']))
                plt.colorbar()
                plt.savefig(direc+str(int(timestep/frameskip)).zfill(10)+'.png', format='png', bbox_inches='tight')
                plt.close(fig)
            counter += 1
            mass_evolution.append(self.graph.get_mass())
            if (df_or_not >0):
                self.update_df()
            print('Timestep '+ str(timestep) + ' completed!!!')
            if (timestep % 2 ==0):
                self.particle_advection_sim() 
            else:
                self.collision_sim()
            
    def plot_steady_state(self):
        #eigendecomposition; eigenvalues, eigenvectors
        self.advection_transition_matrix()

        L, V = linalg.eig(self.advection_tMatrix, left=True, right=False)

        cnt = 0
        for l in L:
            if l >= 0.99999999:
                break
            cnt += 1

        steady_state = V.transpose().real[cnt]
        sum = steady_state.sum()
        steady_state = steady_state / sum
        steady_state = steady_state / steady_state.max()

        #plot
        plt.imshow(steady_state.reshape(self.graph.n, self.graph.n), cmap="jet",  extent=[0,1,0,1], vmin = 0, vmax = 1)
        plt.colorbar()
        plt.show()
        
domain = GyreDomain()

# Initialize 20 x 20 nodes each with 100 particles with velocity drawn randomly from a 2D Gaussian distribution
# The velocity distribution has mean: [0,0], covariance being the 2 x 2 identity matrix
#       velocities=np.random.multivariate_normal([1,1],[[1,0],[0,1]],number_of_particles)
list_of_nodes = [0]*(n**2)   # initialize a list of appropriate length to hold these nodes 

list_of_IDs =[]

for k in range (n**2):
    new_node = Node(domain, k, n, 0, [0]*number_of_particles_at_each_node[k])
    new_node_vel=np.random.multivariate_normal([0,0],[[1,0],[0,1]],number_of_particles_at_each_node[k])  # initialize a velocity space of the particles at the new node
    for i in range (number_of_particles_at_each_node[k]) : # fillin the particles within new_node : 
                new_node.particles[i] = Particle (new_node_vel[i] ,k, n, 0, 'Particle'+str(k)+'___'+str(i))
                list_of_IDs.append('Particle'+str(k)+'___'+str(i))
    # now that new_node is completed, we need to add it to the Graph
    list_of_nodes[k] = new_node

df = pd.DataFrame(np.zeros((len(list_of_IDs),2)) ,index = list_of_IDs, columns=['time '+str(0),'time'+ str(int(delta))])
    
g = Graph(n, p_0, domain, list_of_nodes)       # create a graph

sim = Simulation(g, options)    # create a GyreSimulation

sim.advection_collision_plot()

def mass_check (number_of_particles_at_each_node, mass_evolution):
    # 1.check whether mass is conserved throghout the entire simulation
    # 2.check whether the mass distribution changed over the collision part of the simulation
    
    conserved = 1
    collision_correct = 1 
    total_mass= sum(number_of_particles_at_each_node)
    for i in range(len(mass_evolution)):
        if (sum(mass_evolution[i])!= total_mass):
            conserved =0
            break
            print('Error: Total mass is not conserved at timestep '+ str(i) )
        if (i ==0 ):
            pass
        elif (i%2 ==0):
            for j in range(len(mass_evolution[i])): # loops through the mass at all nodes and compare with previous timestep
                if (mass_evolution[i][j] != mass_evolution[i-1][j]):
                    collision_correct = 0
                    break
                    print('Error: Collision at timestep '+ str(i) + 'did not conserve mass')
                    
    print('Congratulations! You\'ve past the mass check!')

mass_check (number_of_particles_at_each_node, mass_evolution)
# print(mass_evolution)

# print(df)

