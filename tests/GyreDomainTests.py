import unittest

from sig.GyreDomain import *

class TestGyreDomain(unittest.TestCase):
    def test_domain_construction(self):
        domain = GyreDomain()

        # test components of x and y limits
        self.assertAlmostEqual(domain.x_lim[0], 0.0, places=12)
        self.assertAlmostEqual(domain.x_lim[1], 1.0, places=12)
        self.assertAlmostEqual(domain.y_lim[0], 0.0, places=12)
        self.assertAlmostEqual(domain.y_lim[1], 1.0, places=12)

        # test the domain area
        self.assertAlmostEqual(domain.area(), 1.0, places=12)
        self.assertAlmostEqual(domain.initial_mass_density(0.1, 0.25), 1.0, places=12)

    def test_external_velocity(self):
        k = 3 # the number of gyres
        domain = GyreDomain(k, (-1.0, 1.0), (0.0, 2.0))

        # test components of x and y limits
        self.assertAlmostEqual(domain.x_lim[0], -1.0, places=12)
        self.assertAlmostEqual(domain.x_lim[1], 1.0, places=12)
        self.assertAlmostEqual(domain.y_lim[0], 0.0, places=12)
        self.assertAlmostEqual(domain.y_lim[1], 2.0, places=12)

        # test the domain area
        self.assertAlmostEqual(domain.area(), 4.0, places=12)
        self.assertAlmostEqual(domain.initial_mass_density(0.1, 0.1), 0.25, places=12)

        # choose a point in space
        x = 0.34
        y = 0.813

        # evaluate the external velocity field
        vel = domain.external_velocity(x, y)
        self.assertEqual(len(vel), 2)
        self.assertAlmostEqual(vel[0], np.sin(k * np.pi * x) * np.cos(k * np.pi * y), places=12)
        self.assertAlmostEqual(vel[1], - np.cos(k * np.pi * x) * np.sin(k * np.pi *y), places=12)
