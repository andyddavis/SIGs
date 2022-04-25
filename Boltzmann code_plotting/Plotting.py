from Graph import *
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
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18


def PlotGamma(t, graph, direc):
    fig = plt.figure(figsize = (40,30), clear = True)
    Gamma = graph.Gamma()
    ax = plt.gca()
    pc = ax.pcolormesh(graph.x, graph.y, Gamma.T, shading='auto', cmap='jet')
    fig.colorbar(pc)
    plt.savefig(direc+'gamma/'+'fig_gamma-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
    fig.clear()
    plt.close(fig)

def PlotPETE(t, graph, direc):
    # plot PETE
    fig = plt.figure(figsize = (40,30), clear = True)
    PETE = graph.PE_TE()
    ax = plt.gca()
    pc = ax.pcolormesh(graph.x, graph.y, PETE.T, shading='auto', cmap='jet')
    fig.colorbar(pc)
    plt.savefig(direc+'PETE/'+'fig_PETE-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
    fig.clear()
    plt.close(fig)

def PlotKETE(t, graph, direc):
    # plot KETE 
    fig = plt.figure(figsize = (40,30), clear = True)
    KETE = graph.KE_TE()
    ax = plt.gca()
    pc = ax.pcolormesh(graph.x, graph.y, KETE.T, shading='auto', cmap='jet')
    fig.colorbar(pc)
    plt.savefig(direc+'KETE/'+'fig_KETE-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
    fig.clear()
    plt.close(fig)
    
def PlotPE(t, graph, direc):
    # plot PE
    fig = plt.figure(figsize = (40,30), clear = True)
    PE = graph.PE_data()
    ax = plt.gca()
    pc = ax.pcolormesh(graph.x, graph.y, PE.T, shading='auto', cmap='jet')
    fig.colorbar(pc)
    plt.savefig(direc+'PE/'+'fig_PE-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
    fig.clear()
    plt.close(fig)

def PlotTE(t, graph, direc):
    # plot TE 
    fig = plt.figure(figsize = (40,30), clear = True)
    TE = graph.TE_data()
    ax = plt.gca()
    pc = ax.pcolormesh(graph.x, graph.y, TE.T, shading='auto', cmap='jet')
    fig.colorbar(pc)
    plt.savefig(direc+'TE/'+'fig_TE-'+str(t).zfill(6)+'.png', format='png', bbox_inches='tight')
    fig.clear()
    plt.close(fig)

def PlotMass(t, X, Y, massDensity, graph, direc):
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
    
