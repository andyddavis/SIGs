import unittest
import random

from Graph import *

class TestGraph(unittest.TestCase):
    def test_graph_construction(self):
        n = 100 # number of nodes per side
        p0 = 0.1 # parameter that determines the probability of staying at node i

        # by defualt, this uses the zero external velocity field
        domain = Domain()

        graph = Graph(n, p0, domain)
        self.assertEqual(graph.n, n)
        self.assertEqual(len(graph.nodes), n**2)
        self.assertAlmostEqual(graph.p_0, p0, 1.0e-12)

    def test_node_initialisation(self):
        # create graph, empty nodes
        domain = Domain()
        for i in range(0,100):
            n = random.randint(1,100)
            graph = Graph(n, 0.1, domain)
            self.assertEqual(len(graph.nodes), n**2)

            for node, k in zip(graph.nodes, range(n**2)):
                self.assertEqual(node.k, k)
