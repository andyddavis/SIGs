import unittest

from tests.DomainTests import *
from tests.GyreDomainTests import *
from tests.NodeTests import *
from tests.GraphTests import *
from tests.SimulationTests import *
from tests.ParticleTests import *

print("Running unittests...\n")

if __name__=='__main__':
    unittest.main(verbosity=2)
