from Domain import *

class GyreDomain:
    def __init__(self):
        Domain.__init__(self)

    # O---------- define floe movement function (vector field)  ----------O
    def external_velocity(self, x, y):
        return (-y+0.5, x-0.5)
