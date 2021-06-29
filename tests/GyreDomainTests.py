import unittest

from GyreDomain import *

class TestGyreDomain(unittest.TestCase):
    def test_gyre_domain_construction(self):
        domain = GyreDomain()

        self.assertAlmostEqual(domain.x_lim[0], 0.0, 1.0e-12)
        self.assertAlmostEqual(domain.x_lim[1], 1.0, 1.0e-12)
        self.assertAlmostEqual(domain.y_lim[0], 0.0, 1.0e-12)
        self.assertAlmostEqual(domain.y_lim[1], 1.0, 1.0e-12)

    def test_external_velocity(self):
        domain = GyreDomain()

        # choose a point in space
        x = 0.34
        y = 0.813

        # evaluate the external velocity field
        vel = domain.external_velocity(x, y)
        self.assertEqual(len(vel), 2)
        self.assertAlmostEqual(vel[0], -y+0.5, 1.0e-12)
        self.assertAlmostEqual(vel[1], x-0.5, 1.0e-12)
