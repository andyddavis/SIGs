class Domain:
    def __init__(self):
        self.x_lim = (0.0, 1.0)
        self.y_lim = (0.0, 1.0)

    # O---------- define floe movement function (vector field)  ----------O
    def external_velocity(self, x, y):
        return (0.0, 0.0)
