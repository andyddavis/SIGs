from Domain import *

class GyreDomain:
    def __init__(self):
        Domain.__init__(self)

    # define external velocity (vector field)
    def external_velocity(self, x, y):
        return ((y-0.5), -(x-0.5))

    # plot the vector field
    def plot(self):
        x,y = np.meshgrid(np.linspace(self.x_lim[0],self.x_lim[1],11),np.linspace(self.y_lim[0],self.y_lim[1],11))
        (u,v) = self.external_velocity(x,y)
        plt.quiver(x,y,u,v)
        plt.show()
        plt.axis('equal')
