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


nx = 50 # number of nodes in the x direction
ny = 50 # number of nodes in the y direction

Lx = 1.0 # length of the domain in the x direction
Ly = 1.0 # length of the domain in the y direction

m = int(1e5) # the number of particles


def upper_1 (x, y, PE, KE, TE):
    gamma = 0
    if y>=0.5:
        gamma = 1
    return gamma

def lower_1 (x, y, PE, KE, TE):
    gamma = 0
    if y<=0.5:
        gamma = 1
    return gamma

def circle_outside_1(x, y, PE, KE, TE):
    gamma = 0
    if ((x - 0.5)**2 + (y - 0.5)**2) >= 1/9:
        gamma = 1
    return gamma

# parent_direc = '/Users/25016/OneDrive/桌面/0414/KE_TE_jump_discontinuity_0.5/'  #accelaration_0.3_collision_rate_KETE
parent_direc = '/Users/lm4307/Desktop/0414/' 
direc_list = ['upper_1/','lower_1/','circle_outside_1/']
gamma_list = [upper_1,lower_1,circle_outside_1]
steps_list = [5,5,5]
skips_list = [5,5,5] # plot every "skips" time steps
timestep_list = [0.01,0.01,0.01] # the timestep size 

for direc,gamma,nsteps,skips,dt in zip(direc_list,gamma_list,steps_list,skips_list,timestep_list):
    tic()
    direc = parent_direc + direc
    qois = ['gamma', 'PETE', 'KETE', 'mass','PE','TE']
    for qoi in qois:
        if os.path.isdir(direc+qoi):
            shutil.rmtree(direc+qoi)
        os.mkdir(direc+qoi)
    
    # create a regular rectangular graph with m particles
    graph = Graph(nx, Lx, ny, Ly, m, gamma)
    idx = 0
    X =[];Y=[]
    plot_or_not = 0
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
            PlotPE(t, graph, direc)
            PlotTE(t, graph, direc)
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
