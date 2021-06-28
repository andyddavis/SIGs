import statistics as st
import math


class Domain:
    def __init__(self):
        self.x_lim = (0,1)
        self.y_lim = (0,1)

    # O---------- define floe movement function (vector field)  ----------O
    def floe_movement(self,x,y):
        x = x-st.mean(self.x_lim)        # centers the graph at
        y = y-st.mean(self.y_lim)        # the center of limits

        z = (x**2 + y**2)**0.5  # z is magnitude to normalise to unit vectors
        if (z == 0):                # catches division by zero
            z = 1

        return (1,1)           # actual function (unit gyre, clockwise)
