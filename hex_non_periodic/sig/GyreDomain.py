from sig.Domain import *
import numpy as np
class GyreDomain(Domain):
    def __init__(self, direc, numGyres=2, center=(0.0, 0.0), radius = 1.0):
        Domain.__init__(self, direc, center=center, radius=radius)
        self.numGyres = numGyres

    # define external velocity (vector field)
    def external_velocity(self, x, y, time):
        # x_mean = (self.x_lim[0] + self.x_lim[1])/2
        # y_mean = (self.y_lim[0] + self.y_lim[1])/2
        # x_dis = self.x_lim[1] - self.x_lim[0]
        # y_dis = self.y_lim[1] - self.y_lim[0]  
        # return (-x_mean+x, -y_mean+y)
        # return (0*np.sin(self.numGyres * np.pi * x) * np.cos(self.numGyres * np.pi * y),
        #         0*np.cos( self.numGyres* np.pi * x) * np.sin(self.numGyres * np.pi *y))
        # return (np.where(True,np.sin(2*np.pi*time),0), 
        #         np.where(True,np.cos(2*np.pi*time),0))
        return (100,100)
    
        return (np.sin(2*np.pi*x), np.sin(2*np.pi*y))