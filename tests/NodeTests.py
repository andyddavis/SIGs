import unittest

import numpy as np

from sig.Node import *
from sig.GyreDomain import *

class TestNode(unittest.TestCase):
    def test_node_construction(self):
        domain = GyreDomain()

        pos1 = np.array([0.1, 0.1])
        area1 = 0.1
        n1 = Node(domain, pos1, area1)

        self.assertAlmostEqual(n1.position[0], pos1[0], places=12)
        self.assertAlmostEqual(n1.position[1], pos1[1], places=12)

        self.assertAlmostEqual(n1.area, area1, places=12)

        pos2 = np.array([0.25, 0.25])
        area2 = 0.25
        n2 = Node(domain, pos2, area2)

        self.assertAlmostEqual(n2.position[0], pos2[0], places=12)
        self.assertAlmostEqual(n2.position[1], pos2[1], places=12)

        self.assertAlmostEqual(n2.area, area2, places=12)

        # check that the mass density is equal to the initial mass density
        self.assertAlmostEqual(n1.mass_density, domain.initial_mass_density(n1.position[0], n1.position[1])*area1, places=12)
        self.assertAlmostEqual(n2.mass_density, domain.initial_mass_density(n2.position[0], n2.position[1])*area2, places=12)

    def test_add_particle(self):
        domain = GyreDomain()

        pos = np.array([0.1, 0.1])
        area = 0.1
        node = Node(domain, pos, area)

        self.assertEqual(len(node.particles), 0)

        # create a particle to add to the current node
        node.add_particle(np.random.normal(size=2))
        self.assertEqual(len(node.particles), 1)

        # create a particle to add to the current node
        node.add_particle(np.random.normal(size=2))
        self.assertEqual(len(node.particles), 2)

    def test_edge_information(self):
        domain = GyreDomain()

        area = 0.1

        pos = np.array([0.1, 0.1])
        node1 = Node(domain, pos, area)

        pos = np.array([0.15, 0.15])
        node2 = Node(domain, pos, area)

        pos = np.array([0.05, 0.15])
        node3 = Node(domain, pos, area)

        node1.set_neighbors([node1, node2, node3])
        self.assertEqual(len(node1.edges), 2)
        self.assertAlmostEqual(node1.edges[0] [0], domain.directional_distance(node1.position, node2.position) [0], places=12)
        self.assertAlmostEqual(node1.edges[1] [0], domain.directional_distance(node1.position, node3.position) [0], places=12)

    def test_advection_probabilities(self):
        domain = GyreDomain()

        area = 0.1

        pos = np.array([0.1, 0.1])
        node1 = Node(domain, pos, area)

        pos = np.array([0.15, 0.15])
        node2 = Node(domain, pos, area)

        pos = np.array([0.05, 0.15])
        node3 = Node(domain, pos, area)

        pos = np.array([0.05, -0.15])
        node4 = Node(domain, pos, area)

        node1.set_neighbors([node1, node2, node3, node4])
        for i in range(10):
            node1.add_particle(np.random.normal(size=2))

        dt = 0.1
        advection_probabilities = node1.advection_transition_probabilities(dt, advection_scale=2.0)
        self.assertEqual(advection_probabilities.shape[0], len(node1.particles))
        self.assertEqual(advection_probabilities.shape[1], len(node1.neighbors))
        for i in range(len(node1.particles)):
            self.assertAlmostEqual(sum(advection_probabilities[i, :]), 1.0, places=12)
