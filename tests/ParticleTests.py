import unittest
import random

from Particle import *

class TestParticle(unittest.TestCase):
    def test_particle_construction(self):
        # construct a particle
        vel = [1, 1]
        pos = [1, 1]
        mass = 1
        particle = Particle(vel, pos, mass)

        # test particle fields
        self.assertAlmostEqual(particle.vel[0], 1.0, 1.0e-12)
        self.assertAlmostEqual(particle.vel[1], 1.0, 1.0e-12)
        self.assertAlmostEqual(particle.pos[0], 1.0, 1.0e-12)
        self.assertAlmostEqual(particle.pos[1], 1.0, 1.0e-12)
        self.assertAlmostEqual(particle.mass, 1.0, 1.0e-12)
