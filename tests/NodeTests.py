import unittest
import random

from Domain import *
from Node import *

class TestNode(unittest.TestCase):
    def test_node_construction(self):
        domain = Domain()
        n = 5
        k = 1
        node = Node(domain, k, n)

        self.assertEqual(node.n, n)

    def test_node_indices(self):
        domain = Domain()
        # test 100 random nodes to make sure indexing formula works
        for i in range(0,100):
            n = random.randint(1,100)
            k = random.randint(0,n**2 - 1)

            node = Node(domain, k, n)
            self.assertEqual(node.j, k % n, 1.0e-12)
            self.assertAlmostEqual(node.i, (k-node.j)/n, 1.0e-12)

    def test_node_pos(self):
        domain = Domain()
        # rest 100 random nodes to make sure position formula works
        for i in range(0,100):
            n = random.randint(1,100)
            k = random.randint(0,n**2 - 1)
            node = Node(domain, k, n)
            self.assertAlmostEqual(node.pos()[0], (node.j + 0.5) / node.n, 1.0e-12)
            self.assertAlmostEqual(node.pos()[1], 1 - (0.5 + node.i) / node.n, 1.0e-12)

    def test_node_velocity(self):
        domain = Domain()
        # rest 100 random nodes to make sure position formula works
        for i in range(0,100):
            n = random.randint(1,100)
            k = random.randint(0,n**2 - 1)
            node = Node(domain, k, n)
            pos = node.pos()
            self.assertAlmostEqual(node.velocity()[0], domain.external_velocity(pos[0], pos[1])[0], 1.0e-12)
            self.assertAlmostEqual(node.velocity()[1], domain.external_velocity(pos[0], pos[1])[1], 1.0e-12)
