class Particle:                 # Particle class

    def __init__(self, vel, pos, mass):
        self.vel = vel              # velocity (np array)
        self.pos = pos              # position (np array)
        self.mass = mass            # mass

    def __str__(self):
        p = "Position: %.2f, %.2f"% (self.pos[0], self.pos[1])
        v = "\nVelocity: %.2f, %.2f" % (self.vel[0], self.vel[1])
        return(p + v)
