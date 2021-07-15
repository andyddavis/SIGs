#!/usr/bin/env python
# coding: utf-8

# In[87]:


# In this toy model, we represent each ice floe as an object with features: mass, position, velocity.
# For simplicity, we make the following assumptions:
# 1. We assume that after collision, no ice floe will split into smaller pieces.
# 2. We use a 2D model, using only the directions south & north, west & east
# 3. We assume a constant coefficient of friction (COF)
# 4. We ignore wind stress 
# 5. We consider each ice floe as a perfect circle with the same density, in that case, mass of a ice floe is propotional to the square of its radius
# 6. We consider elastic collision 


# In[88]:


COF = 0.2
gravity=9.8
density=1


# In[89]:


import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations


# In[90]:


def distance (x,y):
    return np.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2)


# In[91]:


class Floe:

    def __init__(self, time, radius, number, vel, pos):
        self.time = time
        self.radius = radius
        self.number = number
        self.vel = vel
        self.pos = pos 
        # vel is a numpy array representing velocity
        # pos is a numpy array representing position
        # other attributes are decimal numbers
    def displayinfo (self):
        print('Hello, this is ice floe number %d' % self.number)
        print('Time: %.2f' % self.time)
        print('Mass: %.2f' % self.radius)
        print('Position: %.2f, %.2f' % (self.pos[0], self.pos[1]))
        print('Velocity: %.2f, %.2f' % (self.vel[0], self.vel[1]))
    
    def increase_time(self, delta):
        self.time = self.time + delta
        x=self.pos[0]+delta*self.vel[0]-0.5*gravity*COF*delta*delta
        y=self.pos[1]+delta*self.vel[1]-0.5*gravity*COF*delta*delta
        self.pos=np.array([x,y])
        vx=self.vel[0]-COF*gravity*delta
        vy=self.vel[1]-COF*gravity*delta
        self.vel=np.array([vx,vy])
        


# In[92]:


def collision_or_not (floe1, floe2):
    dist=distance(floe1.pos,floe2.pos)
    size=floe1.radius+ floe2.radius
    return dist<size


# We use the 2D point mass elastic collision formula

# In[93]:


def implement_collision (floe1, floe2):
    m1=floe1.radius**2
    m2=floe2.radius**2
    M=m1+m2
    pos1=np.array(floe1.pos)
    pos2=np.array(floe2.pos)
    vel1=np.array(floe1.vel)
    vel2=np.array(floe2.vel)
    dist_squared=distance(pos1,pos2)**2
    new_vel1=vel1-2*m2/M *np.dot(vel1-vel2, pos1-pos2)/dist_squared * (pos1-pos2)
    new_vel2=vel2-2*m1/M *np.dot(vel2-vel1, pos2-pos1)/dist_squared * (pos2-pos1)
    floe1.vel=new_vel1
    floe2.vel=new_vel2


# In[94]:


def real_collision(floe1,floe2):
    if(collision_or_not (floe1, floe2)):
        implement_collision (floe1, floe2)        


# In[95]:


# Initialize a toy dataset: 

#time, radius, number, vel, pos


F1=Floe(0,1,0,[1,1],[1,1])
F2=Floe(0,2,1,[-2,-1.5],[8,8])
F3=Floe(0,1,2,[-1,1],[1,-2])
F4=Floe(0,0.5,3,[2,2],[-3,-3])
F5=Floe(0,0.1,4,[5,5],[0,0])
set_of_floes=[F1,F2,F3,F4,F5]


# In[96]:


for i in range (100):
    for floe in set_of_floes:
        floe.increase_time(0.1)
    lst=list(itertools.combinations(set_of_floes,r=2))
    for pair in lst:
        real_collision(pair[0],pair[1])
    for i in set_of_floes:
        i.displayinfo()
        print('----------------------------------------')
    print('+++++++++++++++++++++++++++++++++++++++++')
    

