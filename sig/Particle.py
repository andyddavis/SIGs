import numpy as np
import random
import math

class Particle:
    def __init__(self, vel):
        # set the initial velocity
        self.vel = vel

    @staticmethod
    def sample_unit_hypersphere():
        # return a random unit vector
        angle_of_attack = 2*np.pi* random.random()
        return np.array([np.cos(angle_of_attack),np.sin(angle_of_attack)])

    def energy(self):
        return 0.5*np.dot(self.vel, self.vel)

    def inelastic_collision(self, particle2, gamma=1.0):
        pos_diff = Particle.sample_unit_hypersphere()

        # sum of the energies
        nrg = self.energy() + particle2.energy()

        # perfectly elastic scaling
        elastic_scale = -np.dot(pos_diff, self.vel-particle2.vel)

        w = 0.5*(elastic_scale + np.sign(elastic_scale)*np.sqrt(max(0.0, elastic_scale*elastic_scale-4.0*(1.0-gamma)*nrg)))

        self.vel = self.vel + w*pos_diff
        particle2.vel = particle2.vel - w*pos_diff

    def quadratic_drag(self, external_velocity, scale=1.0):
        dum = external_velocity-self.vel
        return scale*np.abs(dum)*dum
