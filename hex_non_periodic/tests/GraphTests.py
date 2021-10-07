import unittest

import numpy as np

from sig.Graph import *
from sig.GyreDomain import *

class TestGraph(unittest.TestCase):
    def test_find_nearest_neighbors(self):
        domain = GyreDomain()

        pos1 = np.array([0.1, 0.1])
        area1 = 0.1
        n1 = Node(domain, pos1, area1)

        pos2 = np.array([0.25, 0.25])
        area2 = 0.25
        n2 = Node(domain, pos2, area2)

        graph = Graph([n1, n2])

        nearest = graph.nearest_node([0.09, 0.11])
        self.assertEqual(nearest, n1)

        nearest = graph.nearest_node([0.29, 0.31])
        self.assertEqual(nearest, n2)

    def test_create_square_lattice(self):
        domain = GyreDomain(x_lim=(-1.0, 1.0), y_lim=(0.0, 2.0))

        n = 10 # the number of nodes on each side
        graph = Graph.create_square_lattice(domain, n)
        self.assertEqual(len(graph.nodes), n*n)

        totalArea = 0.0
        for i in range(n):
            for j in range(n):
                totalArea += graph.nodes[n*i+j].area

                # the position for a [0, 1]x[0, 1] domain
                pos = np.array([(i+0.5)/float(n), (0.5+j)/float(n)])
                # cmpute the position in the domain limits
                pos[0] = (domain.x_lim[1]-domain.x_lim[0])*pos[0]+domain.x_lim[0]
                pos[1] = (domain.y_lim[1]-domain.y_lim[0])*pos[1]+domain.y_lim[0]

                self.assertAlmostEqual(graph.nodes[n*i+j].position[0], pos[0], places=12)
                self.assertAlmostEqual(graph.nodes[n*i+j].position[1], pos[1], places=12)

                # make sure we have stored all of the nearest neighbors
                self.assertEqual(len(graph.nodes[n*i+j].neighbors), 5)
                for neigh, index in zip(graph.nodes[n*i+j].neighbors, [n*i+j, n*((i+1)%n)+j, n*((i-1)%n)+j, n*i+(j+1)%n, n*i+(j-1)%n]):
                    self.assertTrue(index<len(graph.nodes))

                    jneigh = index%n
                    ineigh = (index-jneigh)/n

                    # the position of the neighbor for a [0, 1]x[0, 1] domain
                    pos = np.array([(ineigh+0.5)/float(n), (0.5+jneigh)/float(n)])
                    # cmpute the position in the domain limits
                    pos[0] = (domain.x_lim[1]-domain.x_lim[0])*pos[0]+domain.x_lim[0]
                    pos[1] = (domain.y_lim[1]-domain.y_lim[0])*pos[1]+domain.y_lim[0]

                    self.assertAlmostEqual(neigh.position[0], pos[0], places=12)
                    self.assertAlmostEqual(neigh.position[1], pos[1], places=12)

        self.assertAlmostEqual(totalArea, domain.area(), places=12)
