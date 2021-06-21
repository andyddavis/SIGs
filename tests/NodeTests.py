import unittest

from Domain import *
from Node import *

class TestNode(unittest.TestCase):
    def test_example(self):
        domain = Domain()
        n = 5
        k = 1

        node = Node(domain, k, n)

        self.assertEqual(node.n, n)
