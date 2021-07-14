If you are a Jupyter user, run the .ipynb file by jupyter directly.
If you are using other IDE, please use the py_version file.

Before running the code, there are severl parameters that need to be changed first.

direc is the directory of the folder in which you wish to hold the simulation results.
In my case, '/Users/25016/Desktop/Advection_Collision_Gaussian/', 25016 is the user name of my windows PC.

n is the number of nodes on one side. 
p_0 is the "inertial" parameter (changes probability of staying)
collision_rate is the proportion of particles that will collide within one timestep on a node. Currently, we are using the same collisoin_rate for all nodes at all times.
gamma is the energy ratio for the inelastic collision we are implementing (energy after / energy before), we are also using the same ratio for all collisions
mass_evolution is a 2D array storing the data of the evolution of mass distribution. For each timestep, the mass distribution is stored as a 1D array, which stores the number of particles in node 1, 2, 3, ..., k in ascending order.
plot_or_not: when set to 1, we are plotting the evolution of mass distribution, when set to 0, we are not
df is a dataframe storing the location of all particles at every timestep. It looks like this:
	  time1   time2  time 3
Particle1 3       6      34 
Particle2 32      1      17
Particle3 24      2      99 

number_of_particles_at_each_node is a ist that specifies how many particles there are at each node initially

df_or_not: when set to 1, we are computing df; when set to 0, we are not
Note: Currently, the computation of df is extremely computational expensive!!!

In this simulation, we initialize at particles at n^2 nodes whose velocities follow a 2D Gaussian distribution. Mean: [0,0]  Variance: 2D identity matrix