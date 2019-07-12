""" somewhat flexible plotting library for PCA"""
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#import modules.utils as utils
import scipy as sp
import matplotlib.markers as markers


def plotPCA(V_g, w_g, ctypes, k = 5, fn = None, fig = None, figsize = None,markersize = None):
    if markersize is None:
        markersize =4
    if fig is None:
        if figsize is None:
            fig =  plt.figure(figsize = (k*4, k*4))
        else:
            fig =  plt.figure(figsize = (figsize, figsize))
    ax  = fig.add_subplot(k-1, k-1, 1)

    ### choose colormap and adapt normalization 
    cmap = plt.get_cmap('jet')
    norm = plt.normalize(0, sp.unique(ctypes).shape[0])
    ### plot first k main axes of variation
    for k1 in range(0, k):
        cnt = 1
        for k2 in range(k1 + 1, k):
            ax = fig.add_subplot(k-1, k-1, (k1 * (k-1)) + cnt)
            cnt += 1
            for idx, ct in enumerate(sp.unique(ctypes)):
                c_idx = sp.where(ctypes == ct)[0]
                if c_idx.shape[0] > 0:
                    ax.plot(V_g[k1, c_idx], V_g[k2, c_idx], markers.MarkerStyle.filled_markers[idx % 13], color = cmap(norm(idx)), label = ct, ms = markersize, alpha=0.75)
            ax.set_title('PC %i vs %i' % (k1 + 1, k2 + 1))
            ax.set_xticks(ax.get_xticks()[::2])
            ax.set_yticks(ax.get_yticks()[::2])
            ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)
        if k1 == (k - 1):
            ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.2, 0.5))
    plt.tight_layout()
    if fn is not None:
        plt.savefig(fn+'.pdf', dpi = 1200, format = 'pdf')

