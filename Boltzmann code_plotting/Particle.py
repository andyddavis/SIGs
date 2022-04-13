import numpy as np

# acceleration parameter
acc = 1 # 10,100

class Particle:
    def __init__(self, vel):
        self.vel = vel # the velocity of this particle
        self.eta = 2.0 # the reflection parameter; must be in [1.0, 2.0]

    def Acceleration(self):
        a = np.array([0,0.1])
        a = a-self.vel
        return np.absolute(a)*a*acc

    def UpdateVelocity(self, dt):
        aMean = np.array([0.0]*2)
        aCov = 0.001*np.identity(2)
        self.vel += dt*self.Acceleration() +  np.sqrt(dt)*np.random.multivariate_normal(aMean, aCov)

    def ReflectVelocity(self, normal):
        self.vel -= self.eta*np.dot(normal, self.vel)*normal
