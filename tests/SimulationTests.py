import unittest
import random

from Domain import *
from Graph import *
from Simulation import *

class TestSimulation(unittest.TestCase):

    def test_simulation_construction(self):
        # construct a valid simulation
        domain = Domain()
        t = 100; dt = 0.1
        n = 20; p_0 = dt;
        graph = Graph(n,p_0,domain)
        mass = []
        sim = Simulation(graph,mass,t,dt)
        # test equality/existence of each field
        self.assertEqual(sim.graph, graph)
        self.assertEqual(sim.mass, mass)
        self.assertEqual(sim.time, t)
        self.assertEqual(sim.dt, dt)
        for i in range (0, graph.n**2):
            for j in range (0, graph.n**2):
                self.assertEqual(sim.wind_tMatrix[i,j],0)
                self.assertEqual(sim.idle_tMatrix[i,j],0)

    def test_matrices_initialisation(self):
        # construct a valid simulation
        domain = Domain()
        t = 100; dt = 0.1
        n = 20; p_0 = dt;
        graph = Graph(n,p_0,domain)
        mass = []
        sim = Simulation(graph,mass,t,dt)

        # construct the matrices
        sim.initialise_matrices()
        for node in sim.graph.nodes:
            # test stochasticity
            self.assertAlmostEqual(sum(sim.wind_tMatrix[node.k]), 1.0, 1.0e-12)
            # test p_stay correctness
            self.assertAlmostEqual(sim.wind_tMatrix[node.k,node.k], 1/(1+ p_0 * n * (node.velocity()[0]**2 + node.velocity()[1]**2)**0.5), 1.0e-12)
        pass
