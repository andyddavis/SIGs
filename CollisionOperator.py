import random
import numpy as np

from Particle import *

class CollisionOperator:        # Collision Operator class

    r_1 = 100   # particles (floes) per unit mass
    r_2 = 0.5   # percent of particles colliding
    def __init__(self, mass, v_0):
        self.mass = mass                # mass in the region
        self.v_0 = v_0                  # velocity (u,v) from the node/velocity field
        self.number_of_particles = round(mass * r_1)  # number of particles in region
        self.particles = []             # particles is an numpy array of Particles

    # initialises the particle list # TODO: MAKE THIS FUNCTION OF v_0
    def initialise_particles(self):
        # generate velocities array (2D gaussian), mean (x,y) = v_0, covariance = identity
        velocities=np.random.multivariate_normal([v_0[0], v_0[1]], [[1,0],[0,1]], number_of_particles)
        # generate and add the particles to the particle list
        m_p = mass/number_of_particles          # particles with uniform mass
        for i in range(number_of_particles):
            new_particle = Particle(velocities[i], np.array([0,0]), m_p)
            self.particles.append(new_particle)

    # HELPER METHOD: collides two particles and updates their velocities
    def collide_particles(particle1, particle2):
        # initial velocities & angle of attack
        v1 = particle1.vel                      # particle 1 velocity
        v2 = particle2.vel                      # particle 2 velocity
        theta = 2 * np.pi * random.random()     # angle of attack [0,2pi]

        # ??? # TODO: comment this
        pos_diff = np.array( [np.cos(theta), np.sin(theta)] )   # ????
        value = np.dot(v1 - v2, pos_diff)                       # ???

        # calclate/assign new velocities
        particle1.vel = v1 - value * pos_diff
        particle2.vel = v2 + value * pos_diff

    # collides all particles in pairs and (hence) updates all particle velocities
    def collide(self):
        # find the number of pairs of particles colliding
        number_of_pairs = round(number_of_particles * r_2 / 2)
        index=list(range(len(self.particles)))
        # create the pairs of particles to collide
        pairs=[]
        for i in range(number_of_pairs):
            new_pair = list(random.sample(index, 2))
            pairs.append(new_pair)
            index = list(set(index) - set(new_pair))
        # collide the particles
        for i in pairs:
            collide_particles(self.particles[i[0]], self.particles[i[1]])

    # runs the collision operator and returns the useful values
    def collision_analysis(self):
        K_e = 0     # kinetic energy
        P_e = 0     # potential energy
        v_f = 0     # final mean velocity
        return (K_e, P_e, v_f)      # return kinetic energy, portential energy, and final expected velocity
