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
    def __init__(self, x_lim=(0.0, 1.0), y_lim=(0.0, 1.0)):
        self.x_lim = x_lim
        self.y_lim = y_lim

    def directional_distance(self, point1, point2):
        xlength = self.x_lim[1]-self.x_lim[0]
        ylength = self.y_lim[1]-self.y_lim[0]

        diff = point2-point1
        vecs = [diff, np.array([diff[0]-xlength, diff[1]]), np.array([diff[0], diff[1]-ylength]), np.array([diff[0]-xlength, diff[1]-ylength]), np.array([diff[0]+xlength, diff[1]]), np.array([diff[0], diff[1]+ylength]), np.array([diff[0]+xlength, diff[1]+ylength])]

        vals = [0.0]*len(vecs)
        for i in range(len(vecs)):
            vals[i] = np.linalg.norm(vecs[i])

        ind = vals.index(min(vals))

        return vals[ind], vecs[ind]/vals[ind] if vals[ind]>1.0e-12 else [0.0, 0.0]

    # define external velocity (vector field) function
    def external_velocity(self, x, y):
        return (0.0, 0.0)

    # collision rate function
    def collision_rate(self, x, y):
        return 1.0

    def mean_free_path(self, x, y):
        return 1.0e-2

    def energy_loss_ratio(self, x, y):
        return 0.8

    # domain area
    def area(self):
        return (self.x_lim[1]-self.x_lim[0])*(self.y_lim[1]-self.y_lim[0])

    # initial mass density---defaults to constant
    def initial_mass_density(self, x, y):
        return 1.0/self.area()

    # external velocity plot
    def plot_external_velocity(self, filename='external_velocity', ext='pdf', narrows=15):
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
