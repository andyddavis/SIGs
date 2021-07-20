from sig.GyreDomain import *
from sig.Simulation import *

# create and plot the domain; this defines external functions such as ocean/atmosphere velocity
domain = Domain()
domain.plot_external_velocity()

options = dict()
options["num_nodes"] = 100
options["num_particles"] = 1000000
options["num_timesteps"] = 5000
options["final_time"] = 10.0
sim = Simulation(domain, options)

sim.run()
