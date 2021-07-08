import numpy as np
from scipy import linalg
from scipy import sparse
import math

# import plotting packages and set default figure options
useserif = False # use a serif font with figures?
import matplotlib as mpl
import matplotlib.pyplot as plt
if useserif:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['text.usetex'] = True
else:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams['text.usetex'] = False
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

class Simulation:

    def __init__(self, graph, options):
        self.graph = graph                                              # the graph
        self.options = options                                          # options for the simulation
        self.advection_tMatrix = np.zeros((graph.n**2,graph.n**2))      # alt matrix
        self.diffusion_tMatrix = np.zeros((graph.n**2,graph.n**2))      # diffusion transition matrix (internal forcing)

    # initialises the node advection probabilities and their corresponding indices
    def node_advection_transition_probabilities(self, node):

        # calculate node indices (periodic boundary), ordered: itself, left, right, below, above)
        transition_node_indices = [node.k, (node.j-1) % self.graph.n + self.graph.n * node.i , (node.j+1) % self.graph.n + self.graph.n * node.i, (node.k + self.graph.n ) % self.graph.n**2, (node.k - self.graph.n) % self.graph.n**2]

        # calculate probability of staying
        p_stay = 1 / (1 + self.graph.p_0 * self.graph.n * self.options['timestep length'] * np.linalg.norm(node.velocity()))

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
    def advection_transition_matrix(self):
        for node in self.graph.nodes:
            (indices, probabilities) = self.node_advection_transition_probabilities(node)
            for i in range(0,len(indices)):
                self.advection_tMatrix[node.k, indices[i]] = probabilities[i]

    # initialise the diffusion matrix
    def diffusion_transition_tMatrix(self):
        ## TODO: diffusion matrix
        pass

    def advection_sim(self, mass):
        self.advection_transition_matrix()
        #sparse_tMatrix = sparse.csr_matrix(self.advection_tMatrix)
        counter = 0
        frameskip = 1
        for timestep in range(0, round(self.options['total time']/self.options['timestep length'])):
            # color scheme
            #cmap = 'gray'
            c_map = 'jet'

            # domain
            d = self.graph.domain
            dom = [d.x_lim[0], d.x_lim[1], d.y_lim[0], d.y_lim[1]]

            #mass = mass.dot(sparse_tMatrix)
            # plot
            if(counter % frameskip == 0):
                fig = plt.figure()
                plt.imshow(mass.reshape(self.graph.n,self.graph.n), cmap=c_map, interpolation='nearest',  extent=dom, vmin=0, vmax = 2 * np.sum(mass) / self.graph.n ** 2)
                plt.title("Time: " + str( timestep * self.options['timestep length']))
                plt.colorbar()
                plt.savefig('figures/Step-'+str(int(timestep/frameskip)).zfill(10)+'.png', format='png', bbox_inches='tight')
                plt.close(fig)
            counter += 1
            mass = mass.dot(self.advection_tMatrix)

    def plot_steady_state(self):
        #eigendecomposition; eigenvalues, eigenvectors
        self.advection_transition_matrix()

        L, V = linalg.eig(self.advection_tMatrix, left=True, right=False)

        cnt = 0
        for l in L:
            if l >= 0.99999999:
                break
            cnt += 1

        steady_state = V.transpose().real[cnt]
        sum = steady_state.sum()
        steady_state = steady_state / sum
        steady_state = steady_state / steady_state.max()

        #plot
        plt.imshow(steady_state.reshape(self.graph.n, self.graph.n), cmap="jet",  extent=[0,1,0,1], vmin = 0, vmax = 1)
        plt.colorbar()
        plt.show()
