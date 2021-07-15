#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def initial_check (sim):
    # display info for all particles in the graph 
    for node in sim.graph.nodes:
        for particle in node.particles:
            print('This is particle is at node '+str(particle.k))
            particle.display_info()

