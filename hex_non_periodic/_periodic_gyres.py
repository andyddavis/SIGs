from sig.GyreDomain import *
from sig.Simulation import *
import pandas as pd
import numpy as np
direc = '/Users/25016/OneDrive/桌面/same_weight/'
# create and plot the domain; this efines external functions such as ocean/atmosphere velocity
domain = GyreDomain(direc)

options = dict()
options["num_nodes"] = 20

# num_nodes is the number of nodes on each side

options["num_particles"] = 10000
options["num_timesteps"] =  200
options["final_time"] = 200
options["direc"] = direc
sim = Simulation(domain, options)
sim.run()
# pd.DataFrame(sim.mass_evolution).to_csv(direc+ 'mass.csv')

# pd.DataFrame(sim.PE_evolution).to_csv(direc+ 'PE.csv')

# pd.DataFrame(sim.TE_evolution).to_csv(direc+ 'TE.csv')

# pd.DataFrame(np.array(sim.PE_evolution)/np.array(sim.TE_evolution)).to_csv(direc+ 'ratio.csv')