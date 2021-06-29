import unittest

from Graph import *

class TestGraph(unittest.TestCase):
    def test_graph_construction(self):
        n = 100 # number of nodes per side
        p0 = 0.1 # parameter that determines the probability of staying at node i

        # by defualt, this uses the zero external velocity field
        domain = Domain()

        graph = Graph(n, p0, domain)
        self.assertEqual(graph.n, n)
        self.assertAlmostEqual(graph.p_0, p0, 1.0e-12)
