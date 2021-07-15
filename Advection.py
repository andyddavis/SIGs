import numpy as np
from scipy import linalg


class Advection:

    def __init__(self, graph, dt):
        self.graph = graph
        self.dt = dt
        self.transition_matrix = np.zeros(self.graph.n**2, self.graph.n**2)

    # initialises the node advection probabilities and their corresponding indices
    def node_advection_transition_probabilities(self, node):

        # calculate node indices (periodic boundary), ordered: itself, left, right, below, above)
        transition_node_indices = [node.k, (node.j-1) % self.graph.n + self.graph.n * node.i , (node.j+1) % self.graph.n + self.graph.n * node.i, (node.k + self.graph.n ) % self.graph.n**2, (node.k - self.graph.n) % self.graph.n**2]

        # calculate probability of staying
        p_stay = 1 / (1 + self.graph.p_0 * self.graph.n * dt * np.linalg.norm(node.velocity()))

        # initialise probabilities
        transition_node_probabilities = [0.0] * len(transition_node_indices)

        # if p_stay 1, no need for other calculations (||v|| = 0)
        if abs(p_stay - 1) < 1.0e-13:
            transition_node_probabilities[0] = p_stay
            return transition_node_indices, transition_node_probabilities

        # otherwise, calculate "leaving measure" for each neighbor
        for i in range(1, len(transition_node_indices)):
            # calculate edge vector
            edge = self.graph.nodes[transition_node_indices[i]].pos() - node.pos()
            if (np.linalg.norm(edge) > 0.5):                                           # assuming [0,1] X [0,1] (dependent on x_lim, y_lim - if errors, that's why)
                edge[0] = - np.sign(edge[0]) * (1 / self.graph.n)                      # periodic boundary condition fix (assuming uniform square lattice grid)
                edge[1] = - np.sign(edge[1]) * (1 / self.graph.n)
            # assign  measure
            transition_node_probabilities[i] = max(0.0, np.dot(edge, node.velocity()))

        # normalise probabilies of leaving
        S = sum(transition_node_probabilities)
        transition_node_probabilities = (1 - p_stay) * np.array(transition_node_probabilities) * (1 / S)

        # finally, assign probability of staying and return indices/probabilities
        transition_node_probabilities[0] = p_stay
        return transition_node_indices, transition_node_probabilities

    # create the advection transition matrix
    def initialise_transition_matrix(self):
        for node in self.graph.nodes:
            (indices, probabilities) = self.node_advection_transition_probabilities(node)
            for i in range(0,len(indices)):
                self.transition_matrix[node.k, indices[i]] = probabilities[i]

    # plots the steady state
    def plot_steady_state(self):
        # initialise the advection transition matrix
        self.advection_transition_matrix()
        # eigendecomposition; eigenvalues, eigenvectors
        L, V = linalg.eig(self.transition_matrix, left=True, right=False)

        # find the eigenvalue of 1 (really close to one)
        cnt = 0
        for l in L:
            if l >= 0.99999999:
                break
            cnt += 1

        # normalise the eigenvector corresponding to value of 1
        steady_state = V.transpose().real[cnt]
        sum = steady_state.sum()
        steady_state = steady_state / sum
        steady_state = steady_state / steady_state.max()

        #plot
        plt.imshow(steady_state.reshape(self.graph.n, self.graph.n), cmap="jet",  extent=[0,1,0,1], vmin = 0, vmax = 1)
        plt.colorbar()
        plt.show()
