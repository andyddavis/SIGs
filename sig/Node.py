import numpy as np

from sig.Particle import *

class Node:
    def __init__(self, domain, position, area):
        self.domain = domain
        self.position = position
        self.area = area

        # initially the number of particles is zero
        self.particles = list()
        self.new_particles = list()

        self.mass_density = self.domain.initial_mass_density(self.position[0], self.position[1])*self.area

    def set_neighbors(self, neighbors):
        self.neighbors = neighbors
        # the first neighbor is always itself
        self.compute_edge_information()

    def compute_edge_information(self):
        # the first node is itself
        assert(self.domain.directional_distance(self.position, self.neighbors[0].position) [0]<1.0e-12)

        self.edges = [None]*(len(self.neighbors)-1)
        for i in range(1, len(self.neighbors)):
            self.edges[i-1] = self.domain.directional_distance(self.position, self.neighbors[i].position)

    def add_particle(self, vel):
        self.particles.append(Particle(vel))

    def compute_mass_density(self, total_num_particles):
        return len(self.particles)/(self.area*total_num_particles)

    def compute_expected_velocity(self):
        expected_vel = np.array([0.0, 0.0])
        for particle in self.particles:
            expected_vel += particle.vel
        return expected_vel/(len(self.particles) if len(self.particles)>0 else 1)

    def update_energy(self):
        expected_vel = self.compute_expected_velocity()
        self.kinetic_energy = 0.5*np.dot(expected_vel, expected_vel)

        self.total_energy = 0.0
        for particle in self.particles:
            self.total_energy += 0.5*np.dot(particle.vel, particle.vel)
        self.total_energy /= (len(self.particles) if len(self.particles)>0 else 1)

        self.potential_energy = self.total_energy-self.kinetic_energy

    def acceleration_collision(self, dt):
        # acceleration step
        external_velocity = self.domain.external_velocity(self.position[0], self.position[1])
        for particle in self.particles:
            particle.vel += dt*particle.quadratic_drag(external_velocity)

        # collision step
        gamma = np.exp(0.5*self.mass_density*self.mass_density)/(1.0 + self.mass_density*self.potential_energy)
        epsilon = self.potential_energy/self.total_energy
        #epsilon = self.domain.mean_free_path(self.position[0], self.position[1])
        #gamma = self.domain.energy_loss_ratio(self.position[0], self.position[1])
        collision_rate = self.mass_density*self.domain.collision_rate(self.position[0], self.position[1])/epsilon
        if collision_rate<1.0e-8:
            return

        subtime = 0.0
        while subtime<dt:
            step = min(dt-subtime, 1.0/collision_rate)
            subtime += step

            colliding = random.sample(range(len(self.particles)), int(step*collision_rate*len(self.particles)))
            i = 0
            while i<len(colliding)-1:
                self.particles[colliding[i]].inelastic_collision(self.particles[colliding[i+1]], gamma=gamma)
                i += 2

    def transition_probabilities(self, dt, alpha = 0.5, advection_scale = 1.0, diffusion_scale = 1.0):
        edge_parameter = 0.0
        for edge in self.edges:
            edge_parameter += 1.0/edge[0]
        edge_parameter /= len(self.edges)

        diffusion_transition = alpha*self.diffusion_transition_probabilities(dt, diffusion_parameter = edge_parameter*diffusion_scale)

        advection_transition = (1.0-alpha)*self.advection_transition_probabilities(dt, advection_parameter = edge_parameter*advection_scale)
        for i in range(len(self.particles)):
            advection_transition[i, :] += diffusion_transition

        return advection_transition

    def diffusion_transition_probabilities(self, dt, diffusion_parameter = 1.0):
        prob = np.array([0.0]*len(self.neighbors))

        # compute the probability that we stay
        prob[0] = 1.0/(1.0 + dt*diffusion_parameter*self.mass_density*self.potential_energy)

        # if we leave, there is an equal probability of going to each neighbor
        prob[1:] = np.array([(1.0 - prob[0])/len(self.edges)]*len(self.edges))

        return prob

    def advection_transition_probabilities(self, dt, advection_parameter = 1.0):
        external_velocity = self.domain.external_velocity(self.position[0], self.position[1])
        prob = np.zeros((len(self.particles), len(self.neighbors)))
        for p in range(len(self.particles)):
            # compute the probability that we stay
            #prob[p, 0] = 1.0/(1.0 + dt*advection_parameter*np.linalg.norm(self.particles[p].quadratic_drag(external_velocity)))
            prob[p, 0] = 1.0/(1.0 + dt*advection_parameter*np.linalg.norm(self.particles[p].vel))

            for e in range(len(self.edges)):
                #prob[p, e+1] = max(0.0, np.dot(external_velocity, self.edges[e] [1]) / self.edges[e] [0])
                prob[p, e+1] = max(0.0, np.dot(self.particles[p].vel, self.edges[e] [1]) / self.edges[e] [0])
            prob[p, 1:] = (1.0-prob[p, 0])*prob[p, 1:]/sum(prob[p, 1:])

        return prob
