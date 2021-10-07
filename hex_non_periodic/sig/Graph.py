import numpy as np

from sig.Node import *

class Graph:
    def __init__(self, nodes):
        self.nodes = nodes

    @staticmethod
    def create_hex_lattice(domain, n):
        # an empty array of nodes
        side = domain.radius/n # the side of each small hex
        # we choose the x direction to be pi/3, the y direction to be pi/2
        x_vec = np.array([side*0.5*np.sqrt(3),side*0.5])  
        y_vec = np.array([0.0,side])
        tem = [None]*(2*n+1) # create an empty row
        nodes=[]; neighbors=[]
        for i in range (2*n+1):
            nodes.append(tem.copy())
            neighbors.append(tem.copy())
        # initialize nodes and neighbors to be two empty (2n+1)*(2n+1) list
        # in this graph, we shift the position of each point by:
        # f(x,y) = (x+n, y+n)
            for j in range(2*n+1):
                
                if ( (i-n)>=0 and (j-n)>=0 and ( (i-n+j-n) > n) ): # top right corner is empty
                    continue
                elif ( (i-n)<=0 and (j-n)<=0 and ( (i-n+j-n) < -n)): # bottom left corner is empty
                    continue
                else:
                    pos = np.array([0.0,0.0])+x_vec*(i-n)+y_vec*(j-n)
                    area = domain.area()/(3*float(n)*float(n)+3*float(n)+1) #(there are 3n^2+3n+1 small hexes)
                    nodes[i][j]=Node(domain,pos,area,0)
                    
        # loop through all nodes and store the indexes of its neighbors in the 2d list neighbors
        for i in range(2*n+1):
            for j in range(2*n+1):
                indexes=[[i,j]] # first neighbor is itself
                try:  # neighbor at direction pi/6
                    if(nodes[i+1][j] and (i+1>=0) and (j>=0)):
                        indexes.append([i+1,j])
                except(IndexError):
                    pass
        
                try:  # neighbor at direction pi/2 
                    if(nodes[i][j+1] and (i>=0) and (j+1>=0)):
                        indexes.append([i,j+1])
                except(IndexError):
                    pass
        
                try:  # neighbor at direction 5pi/6
                    if(nodes[i-1][j+1] and (i-1>=0) and (j+1)>=0):
                        indexes.append([i-1,j+1])
                except(IndexError):
                    pass
        
                try:  # neighbor at direction 7pi/6
                    if(nodes[i-1][j] and (i-1>=0) and (j>=0)):
                        indexes.append([i-1,j])
                except(IndexError):
                    pass
        
                try:  # neighbor at direction 3pi/2
                    if(nodes[i][j-1] and (i>=0) and (j-1>=0)):
                        indexes.append([i,j-1])
                except(IndexError):
                    pass
        
                try:  # neighbor at direction 11pi/6
                    if(nodes[i+1][j-1] and (i+1>=0) and (j-1>=0)):
                        indexes.append([i+1,j-1])
                except(IndexError):
                    pass        
                
                neighbors[i][j]=indexes
        # store the neighbors into each node 
        for i in range(2*n+1):
            for j in range(2*n+1):
                if (nodes[i][j]):
                    neighs=[nodes[neigh[0]][neigh[1]]  for neigh in neighbors[i][j]]
                    nodes[i][j].set_neighbors(neighs)                

        return Graph(np.array(nodes))

    # # find the node that is closest to a given point
    # def nearest_node(self, point):
    #     return min(self.nodes, key=lambda node: np.linalg.norm(node.position-point))

    def acceleration_collisions(self, dt):
        for node in self.nodes.flatten():
            if(node):
                node.new_collision(dt)

    # def move_particles(self, dt, initial_PE):
    #     num_particles = 0
    #     for node in self.nodes.flatten():
    #         if (node):
    #             if(len(node.particles)>0):
    #                 num_particles += len(node.particles)
    #                 transition_probs = node.transition_probabilities(initial_PE, dt)
    #                 for p in range(len(node.particles)): # particle_p at "node"
    #                     uni = np.random.uniform()
    #                     cumprob = 0.0
    #                     ind = 0
    #                     while cumprob<uni:
    #                         # try:
    #                         #     #print('Dont have index error')
    #                         #     cumprob += transition_probs[p, ind]
    #                         # except (IndexError):
    #                         #     print('\n Having an index error!\n')
    #                         #     break
    #                         cumprob += transition_probs[p, ind]    
    #                         ind += 1
    #                     node.neighbors[ind-1].new_particles.append(node.particles[p])
    #                 node.particles.clear()


    #     # recompute the mass density and energy at each node
    #     for node in self.nodes.flatten():
    #         if (node):
    #             node.particles, node.new_particles = node.new_particles, node.particles
    #             node.mass_density = node.compute_mass_density(num_particles)
    #             node.update_energy()

    def new_move_particles(self, dt, initial_PE):
        num_particles = 0
        for node in self.nodes.flatten():
            if (node and len(node.particles)>0):
                num_neigh = len(node.neighbors)
                num_particles += len(node.particles)
                transition_probs = node.transition_probabilities(initial_PE, dt)
                for p in range(len(node.particles)):
                    prob= np.cumsum(transition_probs[p,:])  # compute the cumsum of transition prob of particle p
                    cumprob = [0]+list(prob)
                    uni = np.random.uniform()
                    for i in range(num_neigh):
                        # if the random number is between cumprob[i] and cumprob[i+1] 
                        if (uni>=cumprob[i] and uni<=cumprob[i+1]):
                            # move the particle to neighbor i's new_particles
                                node.neighbors[i].new_particles.append(node.particles[p])
                node.particles.clear()  # delete the .particles
        # recompute the mass density and energy at each node
        for node in self.nodes.flatten():
            if (node):
                node.particles, node.new_particles = node.new_particles, node.particles # switch particles with new_particles
                node.mass_density = node.compute_mass_density(num_particles) # recompute mass density
                node.update_energy()



    def get_mass (self):
        # return an array of the mass of flattened nodes
        mass=[]
        for i in self.nodes.flatten():
            if(i):
                mass.append(len(i.particles))
            else:
                mass.append(0)
        return np.array(mass)
    
    def get_PE(self):
        # return an array of the PE of flattened nodes
        PE=[]
        for i in self.nodes.flatten():
            if(i):
                PE.append(i.potential_energy)
            else:
                PE.append(0)
        return np.array(PE)
    
    def get_TE (self):
        # return an array of the total energy of flattened nodes
        TE=[]
        for i in self.nodes.flatten():
            if(i):
                TE.append(i.total_energy)
            else:
                TE.append(0)
        return np.array(TE)
    def increase_time(self, dt):
        # increase the time of all nodes and particles by dt
        for node in self.nodes.flatten():
            if (node):
                node.increase_time(dt)
                for particle in node.particles:
                    particle.increase_time(dt)
                
    def update_energy(self): 
        # recompute energy for each node
        for node in self.nodes.flatten():
            node.update_energy()
 
    def create_scatter_data(self):
        # create data for mass density plot
        # the # of appearance of each node = # of particles at the node        
        X=[]
        Y=[]
        for node in self.nodes.flatten():                                  
            if (node):
                for i in node.particles:
                    X.append(node.position[0])
                    Y.append(node.position[1])
        return (X,Y)
    
    def create_scatter_PE(self):
        # create data for mass density plot
        # the # of appearance of each node = int(PE) at the node        
        X=[]
        Y=[]
        for node in self.nodes.flatten():                                  
            if (node):
                for i in range(int(node.potential_energy)):
                    X.append(node.position[0])
                    Y.append(node.position[1])
        return (X,Y)
        
    def compute_center_position(self):
        # compute position of mass center
        center = np.array([0.0,0.0])
        counter = 0
        for node in self.nodes.flatten():
            if(node):
                for particle in node.particles:
                    center += node.position
                    counter += 1
        return center/counter
    
    def compute_center_PE(self):
        # compute position of mass center (weighted average of the pos of all nodes)
        center = np.array([0.0,0.0])
        total = 0
        for node in self.nodes.flatten():
            if (node):
                PE = node.potential_energy
                center+=PE*node.position
                total+=PE
        return center/total