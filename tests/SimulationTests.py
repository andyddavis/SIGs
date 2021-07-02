import unittest
import random

from Domain import *
from Graph import *
from Simulation import *

class TestSimulation(unittest.TestCase):

    def test_simulation_construction(self):
        # construct a valid simulation
        domain = Domain()
        n = 20; p_0 = 0.1;
        graph = Graph(n,p_0,domain)
        self.assertEqual(len(graph.nodes), n**2)

        dt = 0.1
        options = dict()
        options['timestep length'] = dt

        sim = Simulation(graph, options)

        # test equality/existence of each field
        self.assertEqual(sim.graph, graph)
        self.assertEqual(len(sim.graph.nodes), n**2)
        self.assertEqual(sim.options['timestep length'], dt)
        for i in range (0, graph.n**2):
            for j in range (0, graph.n**2):
                self.assertEqual(sim.advection_tMatrix[i,j],0)

    def test_advection_transition_matrix(self):
        domain = Domain()
        n = 20 # number of nodes per size
        p0 = 0.1 # parameter that determines the probability of staying
        graph = Graph(n, p0, domain)

        options = dict()
        options['timestep length'] = 0.1

        sim = Simulation(graph, options)

        sim.advection_transition_matrix()

    def test_matrices_initialisation(self):
        # construct a valid simulation
        domain = Domain()
        n = 20; p_0 = 0.1;
        graph = Graph(n,p_0,domain)

        options = dict()
        options['timestep length'] = 0.1

        sim = Simulation(graph, options)

        # construct the matrices
        sim.advection_transition_matrix()
        for node in sim.graph.nodes:
            # test stochasticity
            self.assertAlmostEqual(sum(sim.advection_tMatrix[node.k]), 1.0, 1.0e-12)
            # test p_stay correctness
            self.assertAlmostEqual(sim.advection_tMatrix[node.k,node.k], 1/(1+ p_0 * n * (node.velocity()[0]**2 + node.velocity()[1]**2)**0.5), 1.0e-12)
        pass

    def test_node_advection_transition_probabilities(self):
        domain = Domain()
        n = 20; p_0 = 0.1;
        graph = Graph(n,p_0,domain)

        options = dict()
        options['timestep length'] = 0.1

        sim = Simulation(graph, options)

        for node in sim.graph.nodes:
            transition_node_indices, transition_node_probabilities = sim.node_advection_transition_probabilities(node)
            self.assertEqual(len(transition_node_indices),len(transition_node_probabilities))
            for index in range(len(transition_node_probabilities)):
                if (index == 0 ):
                    self.assertAlmostEqual(transition_node_probabilities[index],1,1.0e-12)
                else:
                    self.assertAlmostEqual(transition_node_probabilities[index],0,1.0e-12)
            for index in transition_node_indices:
                neighbor = sim.graph.nodes[index]
                # vertical neighbor check
                if (node.i == 0):
                    self.assertTrue((neighbor.i == 0) or (neighbor.i == 1) or (neighbor.i == n-1))
                if (node.i == n-1):
                    self.assertTrue((neighbor.i == 0) or (neighbor.i == n-2) or (neighbor.i == n-1))
                # interior
                if (node.i != 0 and node.i != n-1):
                    self.assertTrue((node.i == neighbor.i) or (node.i == neighbor.i + 1) or (node.i == neighbor.i-1))
                # horizontal neighbor check
                if (node.j == 0):
                    self.assertTrue((neighbor.j == 0) or (neighbor.j == 1) or (neighbor.j == n-1))
                if (node.j == n-1):
                    self.assertTrue((neighbor.j == 0) or (neighbor.j == n-2) or (neighbor.j == n-1))
                # interior
                if (node.j != 0 and node.j != n-1):
                    self.assertTrue((node.j == neighbor.j) or (node.j == neighbor.j + 1) or (node.j == neighbor.j-1))
