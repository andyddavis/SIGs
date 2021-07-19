import unittest

import numpy as np

from sig.Particle import *

class TestParticle(unittest.TestCase):
    def test_particle_construction(self):
        vel1 = np.array([0.2, 0.0])
        p1 = Particle(vel1)

        vel2 = np.array([1.0, 1.0])
        p2 = Particle(vel2)

        self.assertAlmostEqual(p1.vel[0], vel1[0], 1.0e-12)
        self.assertAlmostEqual(p1.vel[1], vel1[1], 1.0e-12)
        self.assertAlmostEqual(p2.vel[0], vel2[0], 1.0e-12)
        self.assertAlmostEqual(p2.vel[1], vel2[1], 1.0e-12)

    def test_quadratic_drag(self):
        vel = np.array([0.2, 0.0])
        p = Particle(vel)

        external_vel = np.array([0.25, -1.0])

        drag = p.quadratic_drag(external_vel, scale=2.5)

        self.assertAlmostEqual(drag[0], 2.5*abs(external_vel[0]-vel[0])*(external_vel[0]-vel[0]), 1.0e-12)
        self.assertAlmostEqual(drag[1], 2.5*abs(external_vel[1]-vel[1])*(external_vel[1]-vel[1]), 1.0e-12)

    def test_inelastic_collision(self):
        vel1 = np.array([0.2, 0.0])
        p1 = Particle(vel1)

        vec = p1.sample_unit_hypersphere()
        self.assertAlmostEqual(np.linalg.norm(vec), 1.0, places=12)

        self.assertAlmostEqual(p1.energy(), 0.5*np.dot(vel1, vel1), places=12)

        vel2 = np.array([1.0, 1.0])
        p2 = Particle(vel2)

        nrg = p1.energy()+p2.energy()

        # by default, this does elastic collision
        p1.inelastic_collision(p2)
        nrg_after = p1.energy()+p2.energy()
        self.assertAlmostEqual(nrg, nrg_after, places=12)

        # try inelastic collision
        for i in range(10):
            nrg = p1.energy()+p2.energy()

            gamma = random.random()
            p1.inelastic_collision(p2, gamma)
            nrg_after = p1.energy()+p2.energy()
            self.assertTrue(nrg+1.0e-12>=nrg_after)
            self.assertTrue(gamma*nrg-1.0e-12<=nrg_after)
