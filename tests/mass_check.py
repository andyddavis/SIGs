#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def mass_check (number_of_particles_at_each_node, mass_evolution):
    # 1.check whether mass is conserved throghout the entire simulation
    # 2.check whether the mass distribution changed over the collision part of the simulation
    
    conserved = 1
    collision_correct = 1 
    total_mass= sum(number_of_particles_at_each_node)
    for i in range(len(mass_evolution)):
        if (sum(mass_evolution[i])!= total_mass):
            conserved =0
            break
            print('Error: Total mass is not conserved at timestep '+ str(i) )
        if (i ==0 ):
            pass
        elif (i%2 ==0):
            for j in range(len(mass_evolution[i])): # loops through the mass at all nodes and compare with previous timestep
                if (mass_evolution[i][j] != mass_evolution[i-1][j]):
                    collision_correct = 0
                    break
                    print('Error: Collision at timestep '+ str(i) + 'did not conserve mass')
                    
    print('Congratulations! You\'ve past the mass check!')

