import numpy as np

# import plotting packages and set default figure options
useserif = True # use a serif font with figures?
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
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


class Domain:
    """A class that defines all of the functions required to define our sea ice model"""

    def __init__(self, x_lim=(0.0, 1.0), y_lim=(0.0, 1.0)):
        """
        Assumes a square domain with potentially different x and y limits.
        Parameters
        ----------
        x_lim : tuple
            The x limits
        y_lim : tuple
            The y limits
        """
        self.x_lim = x_lim
        self.y_lim = y_lim

    def directional_distance(self, point1, point2):
        """
        Returns the distance between two points and the vector that points between them assuming a doubly-periodic domain.
        Parameters
        ----------
        point1 :
            The first point.
        point2 :
            The second point.

        Returns
        ----------
        The distance between the two points and a unit vector that points from point1 to point2. Moving along that vector is the shortest path from point1 to point2
        """
        xlength = self.x_lim[1]-self.x_lim[0]
        ylength = self.y_lim[1]-self.y_lim[0]

        diff = point2-point1
        vecs = [diff, np.array([diff[0]-xlength, diff[1]]), np.array([diff[0], diff[1]-ylength]), np.array([diff[0]-xlength, diff[1]-ylength]), np.array([diff[0]+xlength, diff[1]]), np.array([diff[0], diff[1]+ylength]), np.array([diff[0]+xlength, diff[1]+ylength])]

        vals = [0.0]*len(vecs)
        for i in range(len(vecs)):
            vals[i] = np.linalg.norm(vecs[i])

        ind = vals.index(min(vals))

        return vals[ind], vecs[ind]/vals[ind] if vals[ind]>1.0e-12 else [0.0, 0.0]

    def external_velocity(self, x, y):
        """The external velocity that the ice feels from the ocean and atmosphere. The default is to return (0,0)."""
        return (0.0, 0.0)

    def collision_rate(self, x, y):
        """The collision rate function that defines the frequently particles collide in the collision operator. Defaults to 1."""
        return 1.0

    def mean_free_path(self, x, y):
        """The nondimensional operator that multiplies the collision operator. Small mean free path means more frequent collisions and larger mean free path means less frequent collisions. Must be between 0 and 1. Defaults to 1e-2."""
        return 1.0e-2

    def energy_loss_ratio(self, x, y):
        """The amount of energy remaining after a two-particle collision. If it is less than one, energy is lost. If it is greater than one, energy is gained. If it is equal to 1, energy is conserved. Defaults to 1."""
        return 1

    def area(self):
        """The area of the entire domain"""
        return (self.x_lim[1]-self.x_lim[0])*(self.y_lim[1]-self.y_lim[0])

    def initial_mass_density(self, x, y):
        """The initial mass density---used to initially distribute the particles. Must integrate to 1"""
        return 1.0/self.area()

    def plot_external_velocity(self, filename='external_velocity', ext='pdf', narrows=15):
        """Plot the external velocity field."""
        x, y = np.meshgrid(np.linspace(self.x_lim[0], self.x_lim[1], narrows), np.linspace(self.y_lim[0], self.y_lim[1], narrows))
        (u, v) = self.external_velocity(x, y)

        fig = plt.figure()
        ax = plt.gca()
        ax.quiver(x, y, u, v)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel(r'$x$')
        ax.set_xlabel(r'$y$')
        plt.savefig(filename+'.'+ext, format=ext, bbox_inches='tight')
        plt.close(fig)
