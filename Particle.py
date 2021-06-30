class Particle:                 # Particle class

    def __init__(self, vel, pos, mass):
        self.vel = vel              # velocity (np array)
        self.pos = pos              # position (np array)
        self.mass = mass            # mass

    def __str__(self):
        p = "Position: %.2f, %.2f\n"% (self.pos[0], self.pos[1])
        v = "Velocity: %.2f, %.2f\n" % (self.vel[0], self.vel[1])
        m = "Mass: %.2f" % (self.mass)
        return(p + v + m)
