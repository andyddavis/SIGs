from sig.Domain import *

class GyreDomain(Domain):
    def __init__(self, numGyres=2, x_lim=(0.0, 1.0), y_lim=(0.0, 1.0)):
        Domain.__init__(self, x_lim=x_lim, y_lim=y_lim)
        self.numGyres = numGyres

    # define external velocity (vector field)
    def external_velocity(self, x, y):
        return (np.sin(self.numGyres * np.pi * x) * np.cos(self.numGyres * np.pi * y), - np.cos( self.numGyres* np.pi * x) * np.sin(self.numGyres * np.pi *y))
