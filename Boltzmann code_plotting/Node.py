import numpy as np
import math

from Particle import *

class Node:
    def __init__(self, x, y, dx, dy, i,j):
        self.x = x # the x location of the node 
        self.y = y # the y location of the node
        self.dx = dx # the length the region in the x direction
        self.dy = dy # the length the region in the y direction
        self.area = self.dx*self.dy # the area of this region
        self.loc = i,j
        # the particles in this region---initially there are none 
        self.particles = list()
        self.newParticles = list()

        self.top = None
        self.bottom = None
        self.left = None
        self.right = None

        self.epsilon = 1.0e-2
        # if self.y<0.5:
        #     self.gamma = 1.0
        # else:
        #     self.gamma = 0.0

    # add a particle in this region
    def CreateParticle(self, vel):
        self.particles.append(Particle(vel))

    def update_energy (self):
        expected_vel = np.array([0.0,0.0])
        TE = 0
        for i in self.particles:
            v = i.vel
            expected_vel += v
            TE += 0.5*np.dot(v,v)
        expected_vel /= len(self.particles) if self.particles else 1
        self.ev = expected_vel
        TE /= len(self.particles) if self.particles else 1
        self.KE = 0.5*np.dot(expected_vel, expected_vel)
        self.TE = TE
        self.PE = TE - self.KE 
        # update gamma 
        gamma = 0
        if self.particles:
            gamma0 = 100
            gamma = (np.tanh(gamma0*(self.PE/self.TE-0.1)) +1)/2
        self.gamma = gamma


    # update the particle positions and velocities
    def ConvectionStep(self, dt):
        if len(self.particles)==0:
            return

        remove = [False]*len(self.particles)
        for p, particle in enumerate(self.particles):
            particle.UpdateVelocity(dt)

            T = 0.0
            while T+1.0e-10<dt and not remove[p]:
                cumprob = np.array([0.0]*4)

                # the probability it goes through the left boundary 
                leftright = self.dy*particle.vel[0]/self.area
                cumprob[0] = max(0.0, -leftright)
                # the probability it goes through the right boundary 
                cumprob[1] = cumprob[0]+max(0.0, leftright)

                # the probability it goes through the top boundary 
                topbot = self.dx*particle.vel[1]/self.area
                cumprob[2] = cumprob[1]+max(0.0, topbot)
                # the probability it goes through the right boundary 
                cumprob[3] = cumprob[2]+max(0.0, -topbot)

                delta = min(dt-T, 1.0/cumprob[3])
                cumprob *= delta
                T += delta

                # the cumulative probability of leaving has to be less than one (decrease the timestep if this fail)
                assert(cumprob[3]-1.0e-10<1.0)

                alpha = np.random.uniform()
                for i in range(4):
                    if alpha<cumprob[i]:
                        if i==0:
                            if self.left==None:
                                particle.ReflectVelocity(np.array([1.0, 0.0]))
                            else:
                                self.left.newParticles += [particle]
                                remove[p] = True
                        elif i==1:
                            if self.right==None:
                                particle.ReflectVelocity(np.array([-1.0, 0.0]))
                            else:
                                self.right.newParticles += [particle]
                                remove[p] = True
                        elif i==2:
                            if self.top==None:
                                particle.ReflectVelocity(np.array([0.0, 1.0]))
                            else:
                                self.top.newParticles += [particle]
                                remove[p] = True
                        elif i==3:
                            if self.bottom==None:
                                particle.ReflectVelocity(np.array([0.0, -1.0]))
                            else:
                                self.bottom.newParticles += [particle]
                                remove[p] = True
                        break

        self.particles = [self.particles[p] for p in range(len(self.particles)) if not remove[p]]

    def CombineParticleLists(self):
        self.particles.extend(self.newParticles)
        self.newParticles = list()

    def CollisionStep(self, mu, dt, idx):
        npart = len(self.particles)
        if npart<2:
            return

        if mu<1.0e-10:
            print(mu, len(self.particles))
        assert(mu>0.0)

        T = 0.0 
        while T+1.0e-10<dt:
            collisionProbability = np.array([0.0]*npart)
            for i in range(npart):
                collisionProbability[i] = self.CollisionRateFuntion()

            maxProb = np.max(collisionProbability)
            assert(maxProb>0.0)

            tau = min(dt-T, self.epsilon/(mu*maxProb))
            collisionProbability *= tau*mu/self.epsilon

            collide = list()
            for i in range(npart):
                if np.random.uniform()<collisionProbability[i]:
                    collide += [i]
            np.random.shuffle(collide)

            for i in range(1, len(collide), 2):
                w = self.SampleUnitHypersphere()
                w *= self.PostCollisionFunction(self.particles[collide[i-1]].vel, self.particles[collide[i]].vel, w, idx)
                self.particles[collide[i-1]].vel += w 
                self.particles[collide[i]].vel -= w

            T += tau
    
    def CollisionRateFuntion(self):
        return 1.0

    def SampleUnitHypersphere(self):
        w = np.random.multivariate_normal(np.array([0.0, 0.0]), np.array([[1.0, 0.0], [0.0, 1.0]]))
        return w/np.linalg.norm(w)

    def PostCollisionFunction(self, v, vprime, w, idx):
        
        we = -w.dot(v-vprime)
        e = 0.5*(np.linalg.norm(v)-np.linalg.norm(vprime))
        return 0.5*(we + math.copysign(1.0, we)*np.sqrt(max(0.0, we*we-4.0*(1.0-self.gamma)*e)))