import numpy as np
import scipy.sparse as sparse
import h5py as h5
import shutil, os
from Node import *
from ttictoc import tic,toc
import time
import pandas as pd
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
m = int(10000) # the number of particles
time = 1000
direc = '/Users/25016/OneDrive/桌面/0414/collision_convergence/10000 particles varying covariance/'  
initialMean = np.array([0.0]*2)
initialCov = np.identity(2)*1
# initialize particles in the nodelis
random = [1.0]
gg = [0.5]
cov = [0.001, 0.01,0.1,1,10,100,1000]
fixed_vel = np.random.multivariate_normal(initialMean, initialCov)
PE_TE_big_list = []; KE_TE_big_list = []; PE_big_list = []; TE_big_list = []; 
for gamma in gg:
    for randomness in random:
        for covariance in cov:
            test_node = Node(0,0,0.02,0.02,0,0)
            PE_TE_list = [];PE_list = []; TE_list = []; KE_TE_list =[];
            for i in range(m):
                if i<randomness*m:
                    test_node.CreateParticle(np.random.multivariate_normal(initialMean, initialCov*covariance))
                else:
                    test_node.CreateParticle(fixed_vel.copy()) 
            for t in range(time):
                test_node.update_energy()
                PE_list.append(test_node.PE)
                TE_list.append(test_node.TE)
                PE_TE_list.append(test_node.PE/test_node.TE)
                KE_TE_list.append(test_node.KE/test_node.TE)
                if t%100 == 0:
                    print('step: ',(t+1),' PE/TE: ',test_node.PE/test_node.TE)
                test_node.NewCollisionStep(gamma)
            PE_TE_big_list.append(PE_TE_list); KE_TE_big_list.append(KE_TE_list)
            PE_big_list.append(PE_list); TE_big_list.append(TE_list)
def plot_and_save(filename, direc, data, plot_type, xlabel, ylabel):
    # plot_type 1: plot
    # plot_type 2: hist
    plt.figure(figsize = (40,30), clear = True)
    plots = []
    for covariance, lst in zip(cov, data):
        if (plot_type == 1):
            plots.append(plt.plot(lst, label = 'cov = '+str(covariance),linewidth = 10))
        if (plot_type == 2):
            plots.append(plt.hist(lst, label = 'cov = '+str(covariance),linewidth = 10))
    plt.xlabel(xlabel,fontsize = 60)
    plt.ylabel(ylabel,fontsize = 60) # ,fontsize = 24
    plt.xticks(fontsize = 50)
    plt.yticks(fontsize = 50)
    plt.title(filename,fontsize = 70)    
    plt.legend(fontsize = 60)     # *plots, names
    plt.savefig(direc + filename)
    plt.close('all')
        
filename1 = 'PE_TE plot with varying covariance'  
filename2 = 'TE_KE plot with varying covariance'  #' (10000 timesteps)'
filename3 = 'PE plot with varying covariance'  
filename4 = 'TE plot with varying covariance'  


step_length = [10,100,1000]

for length in step_length:
    plot_and_save(filename1 + ' (' + str(length) +' timesteps)', direc, 
                  [lst[:length] for lst in PE_TE_big_list], 1, 'time', 'PE/TE')
    plot_and_save(filename2 + ' (' + str(length) +' timesteps)', direc, 
                  [lst[:length] for lst in KE_TE_big_list], 1, 'time', 'KE/TE')
    plot_and_save(filename3 + ' (' + str(length) +' timesteps)', direc, 
                  [lst[:length] for lst in PE_big_list], 1, 'time', 'PE')
    plot_and_save(filename4 + ' (' + str(length) +' timesteps)', direc, 
                  [lst[:length] for lst in TE_big_list], 1, 'time', 'TE')
    # file_name1 = "PE_TE plot" + "_randomness_"+ str(int(100*randomness))
    # file_name2 = "PE_TE hist" + "_randomness_"+ str(int(100*randomness))
    # file_name3 = "PE plot" + "_randomness_"+ str(int(100*randomness))
    # file_name4 = "TE plot" + "_randomness_"+ str(int(100*randomness))
    
    # file_name5 = "(100 timesteps) PE_TE plot" + "_randomness_"+ str(int(100*randomness))
    # file_name6 = "(100 timesteps) PE_TE hist" + "_randomness_"+ str(int(100*randomness))
    # file_name7 = "(100 timesteps) PE plot" + "_randomness_"+ str(int(100*randomness))
    # file_name8 = "(100 timesteps) TE plot" + "_randomness_"+ str(int(100*randomness))
    
    # # plot PE_TE
    # plt.plot(PE_TE_list)
    # plt.xlabel('time')
    # plt.ylabel('PE/TE')
    # plt.savefig(direc + file_name1)
    # plt.close('all')
    
    # # hist PE_TE
    # plt.hist(PE_TE_list)
    # plt.xlabel('PE/TE')
    # plt.ylabel('density')
    # plt.savefig(direc + file_name2)
    # plt.close('all')
    
    # # plot PE
    # plt.plot(PE_list)
    # plt.xlabel('time')
    # plt.ylabel('PE')
    # plt.savefig(direc + file_name3)
    # plt.close('all')
    
    # # plot TE
    # plt.plot(TE_list)
    # plt.xlabel('time')
    # plt.ylabel('TE')
    # plt.savefig(direc + file_name4)
    # plt.close('all')
    
    # # plot PE_TE 100
    # plt.plot(PE_TE_list[:100])
    # plt.xlabel('time')
    # plt.ylabel('PE/TE 100')
    # plt.savefig(direc + file_name5)
    # plt.close('all')
    
    # # hist PE_TE 100
    # plt.hist(PE_TE_list[:100])
    # plt.xlabel('PE/TE 100')
    # plt.ylabel('density')
    # plt.savefig(direc + file_name6)
    # plt.close('all')
    
    # # plot PE 100
    # plt.plot(PE_list[:100])
    # plt.xlabel('time')
    # plt.ylabel('PE 100')
    # plt.savefig(direc + file_name7)
    # plt.close('all')
    
    # # plot TE 100
    # plt.plot(TE_list[:100])
    # plt.xlabel('time')
    # plt.ylabel('TE 100')
    # plt.savefig(direc + file_name8)
    # plt.close('all')
    
    # mean = np.mean(PE_TE_list); variance = np.var(PE_TE_list)
    # data = np.array(  [[mean,variance]+PE_TE_list]) #
    # df = pd.DataFrame(data)
    # df.to_csv(direc+"PE_TE" + "_randomness_"+ str(int(100*randomness))+'.csv')
total_time = toc()
print('Total time for the simulation: ', total_time)
