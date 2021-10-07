import unittest

import numpy as np

from sig.GyreDomain import *
from sig.Simulation import *

class TestSimulation(unittest.TestCase):
    def setUp(self):
        x_lim = (-1.0, 1.0)
        y_lim = (0.0, 2.0)
        domain = GyreDomain(x_lim=x_lim, y_lim=y_lim)

        self.num_particles = 2500

        options = dict()
        options["num_nodes"] = 10
        options["num_particles"] = self.num_particles

        self.sim = Simulation(domain, options)

    def test_simulation_construction(self):
        # check the total number of particles
        cumparts = 0
        for node in self.sim.graph.nodes:
            cumparts += len(node.particles)
        self.assertEqual(cumparts, self.num_particles)

    def test_run_simulation(self):
        self.sim.run()
