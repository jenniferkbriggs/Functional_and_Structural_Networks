# %% 
# This code is made to extract the functional network from Calcium timecourses. 
# Jennifer Briggs 2022

## ---- Options for you to change -------
# set to true if you'd like to see and save figures. set to false if you don't need figures
fig_on = False

# if you want to predefine a savepath. If not, comment out this line by putting at # in front!
savepath = '/Users/briggjen/Documents/GitHub/Functional_and_Structural_Networks'
thr = 0.99


# %% Import packages
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.optimize as op
import easygui #for selecting files using gui
import pandas as pd
import math


# %% Functions and Classes
def lookatnetwork(G):
    degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
    dmax = max(degree_sequence)


    fig = plt.figure("Degree of a random graph", figsize=(8, 8))
    # Create a gridspec for adding subplots of different sizes
    axgrid = fig.add_gridspec(5, 4)

    ax0 = fig.add_subplot(axgrid[0:3, :])
    Gcc = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    pos = nx.spring_layout(Gcc, seed=10396953)
    nx.draw_networkx_nodes(Gcc, pos, ax=ax0, node_size=20)
    nx.draw_networkx_edges(Gcc, pos, ax=ax0, alpha=0.4)
    ax0.set_title("Connected components of G")
    ax0.set_axis_off()

    ax1 = fig.add_subplot(axgrid[3:, :2])
    ax1.plot(degree_sequence, "b-", marker="o")
    ax1.set_title("Degree Rank Plot")
    ax1.set_ylabel("Degree")
    ax1.set_xlabel("Rank")

    ax2 = fig.add_subplot(axgrid[3:, 2:])
    ax2.bar(*np.unique(degree_sequence, return_counts=True))
    ax2.set_title("Degree histogram")
    ax2.set_xlabel("Degree")
    ax2.set_ylabel("# of Nodes")

    fig.tight_layout()

    return degree_sequence #export degree sequence

    plt.savefig(savepath + 'NetworkStats.png')
    plt.clf
    
def loglogcost(x, y):
    # remove any zeros - THIS STEPS ASSUMES THE ONLY THE CONNECTED NODES FIT A POWER LAW
    slope, intercept, r, p, std_err = stats.linregress(np.log10(x), np.log10(y))
    err = (r + 1) #r should be -1 so error is r - (-1)
    return err

def makegraph(thr):
    #Finding optimal threshold based on log log fit
    G = nx.Graph() #make an empty graph
    G.add_nodes_from(list(cor_mat)) #add node for each cell
    #loop over every combination of nodes and add an edge if the correlation is greater than threhsold 
    for x in G.nodes:
        xint = int(x)
        for y in range(xint+1, len(G.nodes)-1): 
            if cor_mat.iloc[xint,y] > thr:
                G.add_edge(x,str(y))

    deg_dict = dict(G.degree()) #export degree of each node into a dictonary
    deg_hist = np.histogram(list(deg_dict.values()))
    deg_hist_bins = list(deg_hist[1])
    deg_hist_vals = list(deg_hist[0])

    x_omit = [i for i,x in enumerate(deg_hist_vals) if x == 0 ]
    
    for i in np.flip(x_omit):
        deg_hist_bins.pop(i)
        deg_hist_vals.pop(i)

    err = loglogcost(deg_hist_bins[1:],deg_hist_vals)

    return err


# %%  Load calcium file
try: # you can directly add the path to your data here
    ca = pd.read_csv('/Users/briggjen/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Anschutz/Islet/TempData/Erli_calcium.csv')
except: #if there is no path, it will ask you to select the folder
    path = easygui.fileopenbox()
    ca = pd.read_csv(path)

try: #if time is in the first axis, we save it and remove
    time = ca.Time
    ca = ca.drop('Time', axis=1)
except:
    print('No time avaliable')
    timeopt = "No"

# %% Compute the correlation matrix
global cor_mat

cor_mat = ca.corr() #computes correlation matrix
if fig_on:
    f = plt.figure(figsize=(19, 15))
    plt.matshow(cor_mat, fignum=f.number)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.title('Correlation Matrix', fontsize=16)
    if 'savepath' not in locals():
        savepath = easygui.fileopenbox()
    plt.savefig(savepath + 'Corrmat.png')
    plt.clf

# %% Computing the network -- need to code in how to find the threshold (8 or power law)

scipy.optimize.minimize(makegraph, 0.9, method = 'Nelder-Mead', bounds = (5, 15) )


if fig_on: 
    from pyvis.network import Network
    net = Network(notebook = True)
    net.from_nx(G)
    net.show(savepath + "network.html")
    deg_seq = lookatnetwork(G)

# # %%
# av_degree = 2*G.number_of_edges()/G.number_of_nodes() #average degree is 2m/n

# # %%
