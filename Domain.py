import numpy as np
import matplotlib.pyplot as plt

class Domain:
    def __init__(self):
        self.x_lim = (0.0, 1.0)
        self.y_lim = (0.0, 1.0)

    # define external velocity (vector field) function
    def external_velocity(self, x, y):
        return (0.0, 0.0)

    # external velocity plot
    def plot(self):
        x,y = np.meshgrid(np.linspace(self.x_lim[0],self.x_lim[1],10),np.linspace(self.y_lim[0],self.y_lim[1],10))
        (u,v) = self.external_velocity(x,y)
        plt.quiver(x,y,u,v)
        plt.show()
        plt.axis('equal')
