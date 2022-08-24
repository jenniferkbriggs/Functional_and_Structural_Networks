# %% 
# This code is made to extract the functional network from Calcium timecourses. 
# Jennifer Briggs 2022

## ---- Options for you to change -------
# set to true if you'd like to see and save figures. set to false if you don't need figures
fig_on = True

# if you want to predefine a savepath. If not, comment out this line by putting at # in front!
global savepath
savepath = '/Users/jkbriggs/Documents/GitHub/Functional_and_Structural_Networks/Examples/'

# How do you want to define the threshold? 
#threshold_opts = 'number_of_connections'
#k = 10
threshold_opts = 'scalefreeish'
min_connect = 5
max_connect = 20


# %% Import packages
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.optimize as op
import easygui #for selecting files using gui
import pandas as pd
import math
import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import askyesno
from Extract_Functional_Net import *


# %%  Load calcium file
try: # you can directly add the path to your data here
    ca = pd.read_csv('/Users/jkbriggs/OneDrive - The University of Colorado Denver/Anschutz/Islet/TempData/Erli_calcium.csv')
except: #if there is no path, it will ask you to select the folder
    path = easygui.fileopenbox('Select Time signal file')
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
        savepath = easygui.diropenbox('Select Folder to Save Data In')
        savepath = savepath + savepath[0] #adds slash
    try: 
        plt.savefig(savepath + 'Corrmat.png')
    except:
        print('Save path is not working')
        savepath = easygui.diropenbox('Select Folder to Save Data In')
        savepath = savepath + savepath[0] #adds slash
        plt.savefig(savepath + 'Corrmat.png')


    plt.clf

# set diagonals equal to zero:
cor_mat = cor_mat.where(cor_mat.values != np.diag(cor_mat),0,cor_mat.where(cor_mat.values != np.flipud(cor_mat).diagonal(0),0,inplace=True))

# %% Computing the network -- need to code in how to find the threshold (8 or power law)

# NOT WORKING
#If how to set threshold is not predefined, choose how to set through gui 
if 'threshold_opts' not in locals():
    root = tk.Tk()

    # click event handler
    def b_degree():
        threshold_opts = 'number_of_connections'
        min_connect = 5
        max_connect = 20

        print('done')
        root.destroy()
        return threshold_opts
    
    def b_scalefree():
        threshold_opts = 'scalefreeish'
        print('done')
        root.destroy()
        return threshold_opts


    top = ttk.Frame(root)
    bottom = ttk.Frame(root)

    top.pack(side=tk.TOP)
    bottom.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

    # create the widgets for the top part of the GUI,
    # and lay them out
    b = ttk.Button(root, text="Predefined average degre", command=b_degree)
    c = ttk.Button(root, text="Scale Free (ish)",  command=b_scalefree)
    b.pack(in_=top, side=tk.LEFT)
    c.pack(in_=top, side=tk.LEFT)

    # start the app
    root.mainloop()


# %%
if threshold_opts == 'number_of_connections':
    # Speficy average number of connections:
    thr = thr_based_on_degree(cor_mat, k)
elif 'scalefreeish':
    maxbnds = float(thr_based_on_degree(cor_mat, min_connect)) #because the minimum connection gives the largest threshold
    minbnds = float(thr_based_on_degree(cor_mat, max_connect))
    bnds = (minbnds, maxbnds)
    x0 = np.mean(bnds)
    # Run optimization
    final_symp = op.minimize(makegraph_err,x0, method = 'Nelder-Mead', bounds = ((minbnds, maxbnds),))
    #find treshold
    thr = final_symp.x[0]

G = makegraph(thr)


if fig_on: 
    from pyvis.network import Network
    net = Network(notebook = True)
    net.from_nx(G)
    net.show(savepath + "network.html")
    deg_seq = lookatnetwork(G)


#save adjacency list
nx.write_adjlist(G, savepath + "Network.adjlist")
# # %%
av_degree = 2*G.number_of_edges()/G.number_of_nodes() #average degree is 2m/n
net_stats = pd.DataFrame
# %%
