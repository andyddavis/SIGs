import numpy as np

from Node import *

class Graph:
    def __init__(self, nx, Lx, ny, Ly, m):
        self.nx = nx # number of nodes in the x direction
        self.ny = ny # number of nodes in the y direction
        self.n = self.ny*self.nx # the total number of nodes

        self.m = m # the number of particles

        self.Lx = Lx # length of the domain in the x direction
        self.Ly = Ly # length of the domain in the y direction

        dx = self.Lx/self.nx # the length of each region in the x direction
        dy = self.Ly/self.ny # the length of each region in the y direction

        # x and y coordinates of the nodes
        self.x = np.linspace(dx/2.0, self.Lx-dx/2.0, self.nx)
        self.y = np.linspace(dy/2.0, self.Ly-dy/2.0, self.ny)

        # create the nodes 
        self.nodes = [None]*self.n
        for i in range(self.n):
            inds = self.LocalIndex(i)
            self.nodes[i] = Node(self.x[inds[0]], self.y[inds[1]], dx, dy, inds[0], inds[1])

        # tell each node who it is connected to
        for i in range(self.nx):
            for j in range(self.ny):
                ind = self.GlobalIndex(i, j)
                if i>0: 
                    self.nodes[ind].left = self.nodes[self.GlobalIndex(i-1, j)]
                if i<self.nx-1:
                    self.nodes[ind].right = self.nodes[self.GlobalIndex(i+1, j)]
                if j>0: 
                    self.nodes[ind].bottom = self.nodes[self.GlobalIndex(i, j-1)]
                if j<self.ny-1:
                    self.nodes[ind].top = self.nodes[self.GlobalIndex(i, j+1)]

        # add particles to the nodes 
        randomness = 0.20 # the randomness of initialization
        initialMean = np.array([0.0]*2)
        initialCov = 0.001*np.identity(2)
        # initialize a list of number of particles at the nodes
        number_of_particles = [0]*self.n
        for s in range(self.m):
            number_of_particles[np.random.randint(0, self.n)] += 1
        for i, number in enumerate(number_of_particles): # at node i, "number" number of particles
            vel = np.random.multivariate_normal(initialMean, initialCov)
            random_number = int(number * randomness) 
            for particle in range(number - random_number):
                self.nodes[i].CreateParticle(vel)
            for random_particle in range(random_number):
                self.nodes[i].CreateParticle(np.random.multivariate_normal(initialMean, initialCov))
        # for i in range(self.m):
        #     self.nodes[np.random.randint(0, self.n)].CreateParticle(np.random.multivariate_normal(initialMean, initialCov))
  
    
  # the global index for each node in the domain
    def GlobalIndex(self, i, j):
        return j*self.nx + i

    # the global index for each node in the domain
    def LocalIndex(self, ind):
        i = ind%self.nx 
        return int(i), int((ind-i)/self.nx)

    def CheckParticleCount(self):
        num = 0 
        for node in self.nodes:
            num += len(node.particles)
        return num==self.m

    def MassDensity(self):
        massDensity = np.array([0.0]*self.n)
        for i in range(self.n):
            massDensity[i] = len(self.nodes[i].particles)/float(self.m)/self.nodes[i].area

        return massDensity.reshape(self.ny, self.nx).T
    
    def Gamma(self):
        Gamma = []
        for node in self.nodes:
            Gamma.append(node.gamma) 

        return np.array(Gamma).reshape(self.ny, self.nx).T
    
    def PE_data(self):
        PE_data = []
        for i in self.nodes:
            PE_data.append(i.PE)

        return np.array(PE_data).reshape(self.ny, self.nx).T

    def TE_data(self):
        TE_data = []
        for i in self.nodes:
            TE_data.append(i.TE)

        return np.array(TE_data).reshape(self.ny, self.nx).T
    
    def KE_data(self):
        KE_data = []
        for i in self.nodes:
            KE_data.append(i.KE)

        return np.array(KE_data).reshape(self.ny, self.nx).T
    
    def PE_TE(self):
        data = []
        for i in self.nodes:
            if i.particles:
                data.append(i.PE/i.TE)
            else:
                data.append(0)

        return np.array(data).reshape(self.ny, self.nx).T
    
    def KE_TE(self):
        data = []
        for i in self.nodes:
            if i.particles:
                data.append(i.KE/i.TE)
            else:
                data.append(0)

        return np.array(data).reshape(self.ny, self.nx).T
    
    # update the particle positions and velocities at each node 
    def ConvectionStep(self, dt):
        # loop through each particle and update the postion/velocity 
        for node in self.nodes:
            node.ConvectionStep(dt)

        for node in self.nodes:
            node.CombineParticleLists()

    def CollisionStep(self, massDensity, dt, idx):
        # loop through each particle and have the particles (potentially) collide
        for i in range(len(self.nodes)):
            inds = self.LocalIndex(i)
            self.nodes[i].CollisionStep(massDensity[inds[0], inds[1]], dt, idx)
