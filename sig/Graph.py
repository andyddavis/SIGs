import numpy as np

from sig.Node import *

class Graph:
    def __init__(self, nodes):
        self.nodes = nodes

    @staticmethod
    def create_square_lattice(domain, n):
        # an empty array of nodes
        nodes = [None]*(n*n)

        neighbors = [None]*(n*n)

        for i in range(n):
            for j in range(n):
                # the position for a [0, 1]x[0, 1] domain
                pos = np.array([(i+0.5)/float(n), (0.5+j)/float(n)])
                # cmpute the position in the domain limits
                pos[0] = (domain.x_lim[1]-domain.x_lim[0])*pos[0]+domain.x_lim[0]
                pos[1] = (domain.y_lim[1]-domain.y_lim[0])*pos[1]+domain.y_lim[0]

                # compute the area associated with this node
                area = domain.area()/(float(n)*float(n))

                # add a new node to the list
                nodes[n*i+j] = Node(domain, pos, area)
                neighbors[n*i+j] = (n*i+j, n*((i+1)%n)+j, n*((i-1)%n)+j, n*i+(j+1)%n, n*i+(j-1)%n)

        for (node, neigh) in zip(nodes, neighbors):
            neighs = [nodes[n] for n in neigh]
            node.set_neighbors(neighs)

        return Graph(nodes)

    # find the node that is closest to a given point
    def nearest_node(self, point):
        return min(self.nodes, key=lambda node: np.linalg.norm(node.position-point))

    def acceleration_collision(self, dt):
        for node in self.nodes:
            node.acceleration_collision(dt)

    def move_particles(self, dt):
        num_particles = 0
        for node in self.nodes:
            num_particles += len(node.particles)
            transition_probs = node.transition_probabilities(dt)

            for p in range(len(node.particles)):
                uni = np.random.uniform()
                cumprob = 0.0
                ind = 0
                while cumprob<uni:
                    cumprob += transition_probs[p, ind]
                    ind += 1
                node.neighbors[ind-1].new_particles.append(node.particles[p])
            node.particles.clear()

        # re compute the mass density and energy at each node
        for node in self.nodes:
            node.particles, node.new_particles = node.new_particles, node.particles
            node.mass_density = node.compute_mass_density(num_particles)
            node.update_energy()
