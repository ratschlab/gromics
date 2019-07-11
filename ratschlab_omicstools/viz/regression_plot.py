"""
Plotting lib for regression prediction
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pdb
import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')

def regression_surface_plot(Y, X1, X2, X3, ols_object, Y_label="Response Var", X1_label="X1", X2_label="X2", file_name="mesh_plot"):
    """surface plot of linear model
    INPUTS
    Y   arr response variable
    X1  arr explanatory variable 1
    X2  arr explanatory variable 2
    X3  arr explanatory variable 3
    ols_object  statsmodels OLS
    Y-label str label
    X1_label    str label
    X2_label    str label
    file_name   str file name
    OUTPUT
    saves a figure
    """
    X1 = X1.reshape(X1.shape[0], 1)
    X2 = X2.reshape(X2.shape[0], 1)
    X3 = X3.reshape(X3.shape[0], 1)
    xx1, xx2 = np.meshgrid(np.linspace(X1.min(), X1.max(), 100), np.linspace(X2.min(), X2.max(), 100))
    Z = ols_object.params[0] * xx1 + ols_object.params[1] * xx2 + ols_object.params[3]
    #image bits
    fig = plt.figure(figsize=(12,8))
    ax = Axes3D(fig) 
    #, azim=-115, elev=15)
    surf = ax.plot_surface(xx1, xx2, Z, cmap=plt.cm.RdBu_r, alpha=0.6, linewidth=0)
    #residuals
    resid = Y - ols_object.predict(np.hstack(( X1, X2, X3, np.ones(X1.shape)))).reshape(Y.shape[0], 1)
    #above 0
    ax.scatter(X1[resid>=0], X2[resid>=0], Y[resid>=0], color='black', alpha=1.0, facecolor='white')
    #below 0
    ax.scatter(X1[resid<0], X2[resid<0], Y[resid<0], color='black', alpha=1.0)
    #labels
    ax.set_xlabel(X1_label)
    ax.set_ylabel(X2_label)
    ax.set_zlabel(Y_label)
    #saving figure
    save_name = '/cbio/grlab/home/dkuo/plots/cqtl_mesh/' + file_name + '.png'
    plt.savefig(save_name, dpi=200)
    print("Saved a mesh plot to %s" % (save_name))


def regression_scatter(Y, X1, X2, file_name="scatter"):
    """scatter plot of linear model
    INPUTS
    Y   arr response variable
    X1  arr explanatary variable 1
    X2  arr explanatary variable 2
    Y_label str label
    X1_label    str label
    X2_label    str label
    file_name   str file suffix
    """
    X1 = X1.reshape(X1.shape[0], 1)
    X2 = X2.reshape(X2.shape[0], 1)
    fig, ax = plt.subplots()
    cax = ax.scatter(X1, X2, c=Y, cmap=plt.cm.YlGnBu, alpha=0.8, lw=0.25, edgecolor='gray')
    cb = plt.colorbar(cax)
    cb.set_label('Normalized Expression')
    ax.set_ylabel('Normalized Methylation')
    ax.set_xlabel('Mutation Status')
    save_name = '/cbio/grlab/home/dkuo/plots/cqtl_mesh/' + file_name + '_scatter.png'
    plt.savefig(save_name, dpi=200)
    print("Saved scatter to %s" % (save_name))
    plt.close()

def regression_scattermat(df, subset_hue = None, file_name="scatter_mat_temp"):
    """scatter matrix of a pandas dataframe
    df          data frame   contains all values
    subset_hue  str          hue for coloring      
    """
    twelve_class = ["#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"]
    sns.set_style('whitegrid')
    if subset_hue in df.columns:
        pair_number = np.unique(df[str(subset_hue)]).shape[0]
        sns.set_palette("Paired", pair_number)
        sns.pairplot(df, hue=subset_hue, size=2.5)
    else:
        sns.pairplot(df, size=2.5)
    save_name = '/cbio/grlab/home/dkuo/plots/cqtl_mesh/' + file_name + '_scatterMatrix.png'
    plt.savefig(save_name, dpi=200)
    print("Saved scatter matrix to %s" % (save_name))
    plt.close()

