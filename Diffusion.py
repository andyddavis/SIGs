import numpy as np

class Diffusion:

    def __init__(self, graph, dt, mass):
        self.graph = graph
        self.transition_matrix = np.zeros(self.graph.n**2, self.graph.n**2)
        self.dt = dt
        self.mass = mass

    def psi(self, k, mass):
        return mass[k] / mass.sum()

    # initialises the node advection probabilities and their corresponding indices
    def node_diffusion_transition_probabilities(self, node, mass):

        # calculate node indices (periodic boundary), ordered: itself, left, right, below, above)
        transition_node_indices = [node.k, (node.j-1) % self.graph.n + self.graph.n * node.i , (node.j+1) % self.graph.n + self.graph.n * node.i, (node.k + self.graph.n ) % self.graph.n**2, (node.k - self.graph.n) % self.graph.n**2]

        # calculate probability of staying
        p_stay = 1 / (1 + self.psi(node.k, mass) * self.options['timestep length'])

        # initialise probabilities
        transition_node_probabilities = [0.0] * len(transition_node_indices)

        # if p_stay 1, no need for other calculations (mass = 0 ?)
        if abs(p_stay - 1) < 1.0e-13:
            transition_node_probabilities[0] = p_stay
            return transition_node_indices, transition_node_probabilities

        # otherwise, calculate "leaving measure" for each neighbor
        for i in range(1, len(transition_node_indices)):
            transition_node_probabilities[i] = (1/4) * self.psi(node.k, mass) * p_stay

        # normalise probabilies of leaving
        S = sum(transition_node_probabilities)
        transition_node_probabilities = (1 - p_stay) * np.array(transition_node_probabilities) * (1 / S)

        # finally, assign probability of staying and return indices/probabilities
        transition_node_probabilities[0] = p_stay
        return transition_node_indices, transition_node_probabilities

    # initialise the diffusion matrix
    def initialise_transition_matrix(self, mass):
        for node in self.graph.nodes:
            (indices, probabilities) = self.node_diffusion_transition_probabilities(node, mass)
            for i in range(0,len(indices)):
                self.diffusion_tMatrix[node.k, indices[i]] = probabilities[i]
