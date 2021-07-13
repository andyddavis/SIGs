class Particle:                 # Particle class

    def __init__(self, vel, mass):
        self.vel = vel              # velocity (np array)
        self.mass = mass            # mass

    def __str__(self):
        v = "Velocity: %.2f, %.2f\n" % (self.vel[0], self.vel[1])
        m = "Mass: %.2f" % (self.mass)
        return(v + m)
