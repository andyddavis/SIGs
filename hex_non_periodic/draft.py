# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 21:42:10 2021

@author: 25016
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sig.Graph import *
from sig.Domain import *
from sig.Node import *
from sig.GyreDomain import *

direc = '/Users/25016/OneDrive/桌面/draft/'
# create and plot the domain; this defines external functions such as ocean/atmosphere velocity
domain = GyreDomain(direc)

side=0.01
x_vec = np.array([side*0.5*np.sqrt(3),side*0.5])
y_vec = np.array([0.0,side])
n=100
tem = [None]*(2*n+1)
nodes=[];  neighbors=[]; nei_count=[]
for i in range (2*n+1):
    nodes.append(tem.copy())
    neighbors.append(tem.copy())
    nei_count.append(tem.copy())
# setup the nodes
for i in range(2*n+1):
    for j in range(2*n+1):
        
        if ( (i-n)>=0 and (j-n)>=0 and ( (i-n+j-n) > n) ):
            continue
        elif ( (i-n)<=0 and (j-n)<=0 and ( (i-n+j-n) < -n)):
            continue
        else:
            pos = np.array([0.0,0.0])+x_vec*(i-n)+y_vec*(j-n)
            area = 3*np.sqrt(3)/2 /(3*float(n)*float(n)+3*float(n)+1)
            nodes[i][j]=Node(domain,pos,area,0)

for i in range(2*n+1):
    for j in range(2*n+1):
        indexes=[[i,j]]
        try:
            if(nodes[i+1][j]):
                indexes.append([i+1,j])
        except(IndexError):
            pass

        try:
            if(nodes[i][j+1] and (i>=0) and (j+1>=0)):
                indexes.append([i,j+1])
        except(IndexError):
            pass

        try:
            if(nodes[i-1][j+1] and (i-1>=0) and (j+1)>=0):
                indexes.append([i-1,j+1])
        except(IndexError):
            pass

        try:
            if(nodes[i-1][j] and (i-1>=0) and (j>=0)):
                indexes.append([i-1,j])
        except(IndexError):
            pass

        try:
            if(nodes[i][j-1] and (i>=0) and (j-1>=0)):
                indexes.append([i,j-1])
        except(IndexError):
            pass

        try:
            if(nodes[i+1][j-1] and (i+1>=0) and (j-1>=0)):
                indexes.append([i+1,j-1])
        except(IndexError):
            pass        
        neighbors[i][j]=indexes


for i in range(2*n+1):
    for j in range(2*n+1):
        if(nodes[i][j]):
            nei_count[i][j]=len(neighbors[i][j])-1
            
# for i in range(2*n+1):
#     for j in range(2*n+1):
#         if (nodes[i][j]):
#             neighs=[ nodes[neigh[0]][neigh[1]]  for neigh in neighbors[i][j]]
#             nodes[i][j].set_neighbors(neighs)                
