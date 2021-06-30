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
                self.assertEqual(sim.wind_tMatrix[i,j],0)
                self.assertEqual(sim.idle_tMatrix[i,j],0)

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
        sim.initialise_matrices()
        for node in sim.graph.nodes:
            # test stochasticity
            self.assertAlmostEqual(sum(sim.wind_tMatrix[node.k]), 1.0, 1.0e-12)
            # test p_stay correctness
            self.assertAlmostEqual(sim.wind_tMatrix[node.k,node.k], 1/(1+ p_0 * n * (node.velocity()[0]**2 + node.velocity()[1]**2)**0.5), 1.0e-12)
        pass
