from Domain import *
import numpy as np

class GyreDomain:
    def __init__(self):
        Domain.__init__(self)

    # define external velocity (vector field)
    def external_velocity(self, x, y):
        x = x
        y = y
        k = 2
        #return (k*y/np.linalg.norm([x,y]), -k*x/np.linalg.norm([x,y]))
        return (np.sin(k * np.pi * x) * np.cos(k * np.pi * y), - np.cos(k * np.pi * x) * np.sin(k * np.pi *y))

    # plot the vector field
    def plot(self):
        x,y = np.meshgrid(np.linspace(self.x_lim[0],self.x_lim[1],11),np.linspace(self.y_lim[0],self.y_lim[1],11))
        (u,v) = self.external_velocity(x,y)
        plt.quiver(x,y,u,v)
        plt.show()
        plt.axis('equal')
