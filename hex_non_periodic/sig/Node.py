import numpy as np
from sig.Particle import *

class Node:
    def __init__(self, domain, position, area, time):
        self.domain = domain
        self.position = position
        self.area = area
        self.particles = list()
        self.new_particles = list()
        self.mass_density = self.domain.initial_mass_density(self.position[0], self.position[1])*self.area
        self.time = time
        
    def __str__(self):
        # print out the location and number of particles of a node
        return 'location: '+ str(self.position[0]) + ', '+ str(self.position[1])+'\nnumber of particles: '+ str(len(self.particles))
    
    def compute_edge_information(self):
        # the first node is itself
        assert(self.domain.directional_distance(self.position, self.neighbors[0].position) [0]<1.0e-12)
        # store the 1/distance for all neighbors in the edges attribute
        self.edges = [None]*(len(self.neighbors)-1)
        for i in range(1, len(self.neighbors)):
            self.edges[i-1] = self.domain.directional_distance(self.position, self.neighbors[i].position)

    def set_neighbors(self, neighbors):
        # stor the neighbors as a list of node as an attribute of node
        self.neighbors = neighbors
        # the first neighbor is always itself
        self.compute_edge_information()

    def compute_mass_density(self, total_num_particles):
        return len(self.particles)/(self.area*total_num_particles)

    def compute_expected_velocity(self):
        
        expected_vel = np.array([0.0, 0.0])
        for particle in self.particles:
            expected_vel += particle.vel
        return expected_vel/(len(self.particles) if len(self.particles)>0 else 1)

    def update_energy(self):
        expected_vel = self.compute_expected_velocity()
        self.kinetic_energy = 0.5*np.dot(expected_vel, expected_vel)

        self.total_energy = 0.0
        for particle in self.particles:
            self.total_energy += 0.5*np.dot(particle.vel, particle.vel)
        self.total_energy /= (len(self.particles) if len(self.particles)>0 else 1)

        self.potential_energy = self.total_energy-self.kinetic_energy

    def new_collision(self, dt):
        external_velocity = self.domain.external_velocity(self.position[0], self.position[1], self.time)
        for particle in self.particles:
            particle.vel += dt * particle.quadratic_drag(external_velocity)
            
        collision_rate = self.potential_energy / self.total_energy
        if (collision_rate < 1.0e-3):
            return
        index=list(range(len(self.particles)))
        number_of_pairs=int(len(self.particles)*collision_rate/2)
        pairs=list(range(number_of_pairs)) # initialize pairs to be of equal length
        # Record the pairs of indexes of particles to collide in the list pairs
        for i in range(number_of_pairs):
            new_pair=list(random.sample(index,2))
            pairs[i]=new_pair
            index=list(set(index)-set(new_pair))
        for i in pairs:
            self.particles[i[0]].inelastic_collision(self.particles[i[1]])


    def acceleration_collision(self, dt):
        # acceleration step
        external_velocity = self.domain.external_velocity(self.position[0], self.position[1],self.time)
        for particle in self.particles:
            particle.vel += dt*particle.quadratic_drag(external_velocity)

        # collision step
        self.update_energy()
        epsilon = self.domain.mean_free_path(self.position[0], self.position[1])
        gamma = self.domain.energy_loss_ratio(self.position[0], self.position[1])
        if (len(self.particles)<2):
            return
        collision_rate = self.mass_density * (self.potential_energy / self.total_energy) /epsilon
        # self.domain.collision_rate(self.position[0], self.position[1])
        if (collision_rate < 1.0e-3):
            return
        subtime = 0.0
        while subtime<dt:
            step = min(dt-subtime, 1.0/collision_rate)
            subtime += step

            colliding = random.sample(range(len(self.particles)), int(step*collision_rate*len(self.particles)))
            i = 0
            while i<len(colliding)-1:
                self.particles[colliding[i]].inelastic_collision(self.particles[colliding[i+1]], gamma=gamma)
                i += 2

    def transition_probabilities(self,initial_PE, dt, alpha = 1.0, advection_scale = 1.0, diffusion_scale = 1.0 ):
        edge_parameter = 0.0
        for edge in self.edges:
            edge_parameter += 1.0/edge[0]
        edge_parameter /= len(self.edges)
        diffusion_transition = alpha*self.new_diffusion_transition_probabilities(initial_PE, 
                                        dt ) #diffusion_parameter = edge_parameter*diffusion_scale

        advection_transition = (1.0-alpha)*self.advection_transition_probabilities(dt,
                                    advection_parameter = edge_parameter*advection_scale)
        for i in range(len(self.particles)):
            advection_transition[i, :] += diffusion_transition
        return advection_transition

    # def diffusion_transition_probabilities(self, dt, diffusion_parameter = 1.0):
    #     # compute the probability that we stay
    #     prob = np.array([0.0]*len(self.neighbors))
    #     PE = self.potential_energy*len(self.particles)
    #     list_of_PE =[i.potential_energy*len(i.particles) for i in self.neighbors[1:]]
    #     value=max(0,PE-min(list_of_PE))
    #     list_of_PE.append(PE)
    #     if (  (max(list_of_PE)-min(list_of_PE)) > 1.0e-6 ):
    #         psi=value/(max(list_of_PE)-min(list_of_PE))
    #         p_stay=1/(1+dt*psi)
    
    #         prob[0]=p_stay
            
    #         # assign a weighted probability of leaving 
    #         weight = list()
    #         weight = [max(0,PE-len(neigh.particles)*neigh.potential_energy) for neigh in self.neighbors[1:]]
    #         # normalize the weight
    #         Sum = sum(weight)
    #         if (Sum < 1.0e-8):
    #             return prob
    #         weight=np.array(weight)/Sum
    #         # if we leave, there is an equal probability of going to each neighbor
    #         prob[1:] = weight*(1.0 - prob[0])
    #         #np.array([(1.0 - prob[0])/len(self.edges)]*len(self.edges))
    #         return prob
        
    #     else:
    #         prob[0]=1.0
    #         return prob
        
    def new_diffusion_transition_probabilities(self,initial_PE, dt, diffusion_parameter =1.0 ):
        prob = np.array([0.0]*len(self.neighbors))
        PE = self.potential_energy*len(self.particles)
        list_of_PE = [i.potential_energy*len(i.particles)  for i in self.neighbors[1:]]
        list_of_weight = [max(0.0,PE-PE_j)/initial_PE  for PE_j in list_of_PE]
        maximum=max(list_of_weight)
        list_of_weight=[maximum for i in list_of_weight]
        x = max(list_of_weight)*dt*diffusion_parameter
        p_leave = 2/(1+np.exp(-x))-1
        if (p_leave < 1.0e-8):
            prob[0]=1.0
            assert (abs(1-sum(prob))<1.0e-8)
            return prob
        else:
            prob[0] = 1-p_leave
            weight = np.array(list_of_weight)/sum(list_of_weight)
            prob[1:] = weight*p_leave
            try:
                assert (abs(1-sum(prob))<1.0e-8)
            except:
                print("The total prob is ")
                print(sum(prob))
                quit()
            return prob
    def advection_transition_probabilities(self, dt, advection_parameter = 1.0):
        #external_velocity = self.domain.external_velocity(self.position[0], self.position[1], self.time)
        prob = np.zeros((len(self.particles), len(self.neighbors)))
        for p in range(len(self.particles)):
            # compute the probability that we stay
            #prob[p, 0] = 1.0/(1.0 + dt*advection_parameter*np.linalg.norm(self.particles[p].quadratic_drag(external_velocity)))
            prob[p, 0] = 1.0/(1.0 + dt*advection_parameter*np.linalg.norm(self.particles[p].vel))
            for e in range(len(self.edges)):
                #prob[p, e+1] = max(0.0, np.dot(external_velocity, self.edges[e] [1]) / self.edges[e] [0])
                prob[p, e+1] = max(0.0, np.dot(self.particles[p].vel, self.edges[e] [1]) / self.edges[e] [0])
            prob[p, 1:] = (1.0-prob[p, 0])*prob[p, 1:]/sum(prob[p, 1:])
        return prob

    def increase_time (self, dt):
        self.time = self.time + dt
        
    def add_particle(self, vel,time):
        self.particles.append(Particle(vel, time))