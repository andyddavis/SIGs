import unittest

from sig.Domain import *

class TestDomain(unittest.TestCase):
    def test_domain_construction(self):
        domain = Domain()

        # test components of x and y limits
        self.assertAlmostEqual(domain.x_lim[0], 0.0, places=12)
        self.assertAlmostEqual(domain.x_lim[1], 1.0, places=12)
        self.assertAlmostEqual(domain.y_lim[0], 0.0, places=12)
        self.assertAlmostEqual(domain.y_lim[1], 1.0, places=12)

        # test the domain area
        self.assertAlmostEqual(domain.area(), 1.0, places=12)

        p1 = np.array([0.1, 0.2])
        p2 = np.array([0.1, 0.2])
        dist, vec = domain.directional_distance(p1, p2)
        self.assertAlmostEqual(dist, 0.0, places=12)
        self.assertAlmostEqual(vec[0], 0.0, places=12)
        self.assertAlmostEqual(vec[1], 0.0, places=12)

        p1 = np.array([0.1, 0.2])
        p2 = np.array([0.15, 0.3])
        dist, vec = domain.directional_distance(p1, p2)
        self.assertAlmostEqual(dist, np.linalg.norm(p2-p1), places=12)
        self.assertAlmostEqual(vec[0], (p2[0]-p1[0])/dist, places=12)
        self.assertAlmostEqual(vec[1], (p2[1]-p1[1])/dist, places=12)

        p1 = np.array([0.9, 0.2])
        p2 = np.array([0.05, 0.2])
        dist, vec = domain.directional_distance(p1, p2)
        self.assertAlmostEqual(dist, 0.15, places=12)
        self.assertAlmostEqual(vec[0], 1.0, places=12)
        self.assertAlmostEqual(vec[1], 0.0, places=12)

        dist, vec = domain.directional_distance(p2, p1)
        self.assertAlmostEqual(dist, 0.15, places=12)
        self.assertAlmostEqual(vec[0], -1.0, places=12)
        self.assertAlmostEqual(vec[1], 0.0, places=12)

        p1 = np.array([0.5, 0.95])
        p2 = np.array([0.5, 0.05])
        dist, vec = domain.directional_distance(p1, p2)
        self.assertAlmostEqual(dist, 0.1, places=12)
        self.assertAlmostEqual(vec[0], 0.0, places=12)
        self.assertAlmostEqual(vec[1], 1.0, places=12)

    def test_external_velocity(self):
        domain = Domain()

        # choose a point in space
        x = 0.34
        y = 0.813

        # evaluate the external velocity field
        vel = domain.external_velocity(x, y)
        self.assertEqual(len(vel), 2)
        self.assertAlmostEqual(vel[0], 0.0, places=12)
        self.assertAlmostEqual(vel[1], 0.0, places=12)
