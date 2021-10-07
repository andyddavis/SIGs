#!/usr/bin/env python
# coding: utf-8

# In[7]:


def time_check(sim, t):
# check if the time after simulation is correct
# sim is the simulation object, t is the total time of simulation
# loop through all partilces and check if the time attribute agrees with the previously specified total time

    correct = 1
    for node in sim.graph.nodes:
        for particle in node.particles:
            if (particle.time != t):
                correct=0
                break
    if (correct):
        print('Congratulations: You\'ve past all time check')
    else:
        print('Error: there are some bugs in the time update of the simulation')

