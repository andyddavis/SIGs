import numpy as np
import scipy.sparse as sparse

import h5py as h5

import os

from Graph import *


from ttictoc import tic,toc

import time

tic()
# import plotting packages and set default figure options
useserif = True # use a serif font with figures?
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
if useserif:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['text.usetex'] = False
else:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams['text.usetex'] = False
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

nx = 50 # number of nodes in the x direction
ny = 50 # number of nodes in the y direction

Lx = 1.0 # length of the domain in the x direction
Ly = 1.0 # length of the domain in the y direction

m = int(1e5) # the number of particles

direc = '/Users/25016/OneDrive/桌面/simu0311/copy/'  #accelaration_0.3_collision_rate_KETE
# direc = '/Users/lm4307/Desktop/simu0311/p_1.5/' 

# create a regular rectangular graph with m particles
graph = Graph(nx, Lx, ny, Ly, m)
dt = 0.01 # the timestep size 
idx = 0
nsteps = 1500
X =[];Y=[]
plot_or_not = 0

skips = 5 # plot every "skips" time steps

# get the location of the nodes
for node in graph.nodes:
        X.append(node.x); Y.append(node.y)
for t in range(nsteps):
    if t%skips ==0:
        plot_or_not =1
    else:
        plot_or_not = 0

    print('timestep', t+1, 'of', nsteps)
    for node in graph.nodes:
        node.update_energy()
        
    if plot_or_not:
        
        # plot gamma 
        fig1 = plt.figure(figsize = (40,30), clear = True)
        Gamma = graph.Gamma()
        ax = plt.gca()
        pc = ax.pcolormesh(graph.x, graph.y, Gamma.T, shading='auto', cmap='jet')
        fig1.colorbar(pc)
        plt.savefig(direc+'gamma/'+'fig_gamma-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
        fig1.clear()
        plt.close(fig1)
        # plot PETE
        fig2 = plt.figure(figsize = (40,30), clear = True)
        PETE = graph.PE_TE()
        ax = plt.gca()
        pc = ax.pcolormesh(graph.x, graph.y, PETE.T, shading='auto', cmap='jet')
        fig2.colorbar(pc)
        plt.savefig(direc+'PETE/'+'fig_PETE-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
        fig2.clear()
        plt.close(fig2)
        # plot KETE 
        fig3 = plt.figure(figsize = (40,30), clear = True)
        KETE = graph.KE_TE()
        ax = plt.gca()
        pc = ax.pcolormesh(graph.x, graph.y, KETE.T, shading='auto', cmap='jet')
        fig3.colorbar(pc)
        plt.savefig(direc+'KETE/'+'fig_KETE-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
        fig3.clear()
        plt.close(fig3)
    print('\tstarting convection step.')
    graph.ConvectionStep(dt)
    print('\tfinished convection step.')
    print('\tstarting collision step.')
    massDensity = graph.MassDensity()
    graph.CollisionStep(massDensity, dt, idx)
    print('\tfinished collision step.')
    if plot_or_not:
        # plot mass
        
        fig = plt.figure(figsize = (40,30), clear = True)
        ax = plt.gca()
        pc = ax.pcolormesh(graph.x, graph.y, massDensity.T, shading='auto', cmap='jet')
        # print expected velocities
        U=[];V=[]
        for node in graph.nodes:
            U.append(node.ev[0]); V.append(node.ev[1])
        plt.quiver(X,Y,U,V)
        #pc = ax.pcolormesh(graph.x, graph.y, massDensity.T, shading='gouraud', cmap='jet')
        fig.colorbar(pc)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylabel(r'$y$')
        ax.set_xlabel(r'$x$')
        #ax.set_xlim([x[0], x[-1]])
        #ax.set_ylim([0, 1.1*max(xmarginal)])
        
        plt.savefig(direc+'mass/'+'fig_mass-density-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
        fig.clear()
        plt.close(fig)
        
    plt.close('all')

total_time = toc()
print('Total time for the simulation: ', total_time)
