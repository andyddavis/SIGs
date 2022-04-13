import numpy as np
import scipy.sparse as sparse

import h5py as h5

import shutil, os

from Graph import *
from Plotting import *

from ttictoc import tic,toc

import time

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


tic()

nx = 50 # number of nodes in the x direction
ny = 50 # number of nodes in the y direction

Lx = 1.0 # length of the domain in the x direction
Ly = 1.0 # length of the domain in the y direction

m = int(1e5) # the number of particles

#direc = '/Users/25016/OneDrive/桌面/simu0311/copy/'  #accelaration_0.3_collision_rate_KETE
#direc = '/Users/lm4307/Desktop/simu0311/p_1.5/' 
direc = 'figures/'

qois = ['gamma', 'PETE', 'KETE', 'mass']
for qoi in qois:
    if os.path.isdir('figures/'+qoi):
        shutil.rmtree('figures/'+qoi)
    os.mkdir('figures/'+qoi)

# create a regular rectangular graph with m particles
graph = Graph(nx, Lx, ny, Ly, m)
dt = 0.01 # the timestep size 
idx = 0
nsteps = 1500
X =[];Y=[]
plot_or_not = 0

skips = 1 # plot every "skips" time steps

# get the location of the nodes
for node in graph.nodes:
        X.append(node.x); Y.append(node.y)

# the initial mass density
massDensity = graph.MassDensity()

for t in range(nsteps):
    if t%skips ==0:
        plot_or_not =1
    else:
        plot_or_not = 0

    print('timestep', t+1, 'of', nsteps)
    for node in graph.nodes:
        node.update_energy()
        
    if plot_or_not:        
        PlotGamma(t, graph, direc)
        PlotPETE(t, graph, direc)
        PlotKETE(t, graph, direc)
        PlotMass(t, X, Y, massDensity, graph, direc)
        plt.close('all')

    print('\tstarting convection step.')
    graph.ConvectionStep(dt)
    print('\tfinished convection step.')
    print('\tstarting collision step.')
    # update the mass density
    massDensity = graph.MassDensity()
    graph.CollisionStep(massDensity, dt, idx)
    print('\tfinished collision step.')

total_time = toc()
print('Total time for the simulation: ', total_time)
