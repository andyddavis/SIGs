import unittest

from Domain import *

class TestDomain(unittest.TestCase):
    def test_domain_construction(self):
        domain = Domain()

        self.assertAlmostEqual(domain.x_lim[0], 0.0, 1.0e-12)
        self.assertAlmostEqual(domain.x_lim[1], 1.0, 1.0e-12)
        self.assertAlmostEqual(domain.y_lim[0], 0.0, 1.0e-12)
        self.assertAlmostEqual(domain.y_lim[1], 1.0, 1.0e-12)

    def test_external_velocity(self):
        domain = Domain()

        # choose a point in space
        x = 0.34
        y = 0.813

        # evaluate the external velocity field
        vel = domain.external_velocity(x, y)
        self.assertEqual(len(vel), 2)
        self.assertAlmostEqual(vel[0], 0.0, 1.0e-12)
        self.assertAlmostEqual(vel[1], 0.0, 1.0e-12)