import numpy as np

# import plotting packages and set default figure options
useserif = True # use a serif font with figures?
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
if useserif:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['text.usetex'] = False
else:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams['text.usetex'] = False
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

direc='/Users/25016/Desktop/outward_radial/'  



class Domain:
    def __init__(self, direc, center=(0.0, 0.0), radius = 1.0):
        # radius is the length of the side of the large hex
        
        self.center = center
        self.radius = radius
        self.direc = direc

    def directional_distance(self, point1, point2):
        # compute the distance between two points, return the norm and unit vector
        distance = point2 - point1
        val =np.linalg.norm(distance)
        return val, distance/val if val>1.0e-12 else [0.0, 0.0]

    # define external velocity (vector field) function
    def external_velocity(self, x, y, time):
        return (1.0, 0.0)

    # collision rate functionn
    def collision_rate(self, x, y):
        return (0.5 if x<0.5 else 1.5)

    def mean_free_path(self, x, y):
        return 1.0e-2

    def energy_loss_ratio(self, x, y):
        return 0.8

    # domain area
    def area(self):
        return (3*np.sqrt(3)/2 * self.radius * self.radius)

    # initial mass density---defaults to constant
    def initial_mass_density(self, x, y):
        return 1.0/self.area()

    # external velocity plot
    def plot_external_velocity(self, time, dt, zfill, filename='external_velocity', ext='pdf', narrows=15):
        x, y = np.meshgrid(np.linspace(-self.radius, self.radius, narrows),
                           np.linspace(-self.radius, self.radius, narrows))
        (u, v) = self.external_velocity(x, y, time)

        fig = plt.figure()
        ax = plt.gca()
        ax.quiver(x, y, u, v)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel(r'$x$')
        ax.set_xlabel(r'$y$')
        plt.savefig(self.direc + filename+ str(round(time/dt)).zfill(zfill)+'.'+ext, format=ext, bbox_inches='tight')
        plt.close(fig)
