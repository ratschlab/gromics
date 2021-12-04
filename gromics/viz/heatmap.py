import matplotlib
import numpy as np
import scipy.cluster.hierarchy as spch
import scipy.spatial.distance as spsd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import pdb

def cluster(mat, distance='euclidean', method='single', dim1=True, dim2=True):
    """
    This function takes a matrix and clusters it in the given dimensions.

    Input:
        mat - input matrix 2D
        distance - distance measure passed to scipy.cluster.hierarchy as keyword metric
                   default is: euclidean
        method - clustering method passed to scipy.cluster.hierarchy as keyword method
                 default is: single
        dim1 - bool that tells to cluster along first dimension (default: true)
        dim2 - bool that tells to cluster along second dimension (default: true)

    Output:
        Returns a 4 tuple containing the permuted idx for the two dimensions and the linkage objects.
            lvs1 - permuted idx for dim1
            lv22 - permuted idx for dim2
            lnk1 - linkage for dim1
            lnk2 - linkage for dim2
        The tuple is returned in that order. If the respective dimension was not to be clustered (dim1/2=False), then
        the permuted idx is just the sorted range and the lvs is None.
    """

    if dim1 and mat.shape[0] > 0:
        print('Compute %s clustering in first dimension' % distance)
        if distance == 'nandist':
             dist = _nanDist(mat)
        else:
             dist = spsd.pdist(mat, distance)
             dist = spsd.squareform(dist)
     
        lnk1 = spch.linkage(dist, method=method, metric=distance)
        dendro = spch.dendrogram(lnk1, p=100000, no_plot=True, truncate_mode='mtica')
        lvs = dendro['leaves']
    else:
        lvs = np.arange(mat.shape[0])
        lnk1 = None
 
    if dim2 and mat.shape[0] > 1:
        print('Compute %s clustering in second dimension' % distance)
        if distance == 'nandist':
            dist = _nanDist(mat.T)
        else:
            dist = spsd.pdist(mat.T, distance)
            dist = spsd.squareform(dist)
 
        lnk2 = spch.linkage(dist, method=method, metric=distance)
        dendro = spch.dendrogram(lnk2, p=100000, no_plot=True, truncate_mode='mtica')
        lvs2 = dendro['leaves']
    else:
        lvs2 = np.arange(mat.shape[1])
        lnk2 = None

    return (lvs, lvs2, lnk1, lnk2)


def makeHeatmapCluster(mat, fn = None, tit = None, xlab = None, ylab = None, cmap = cm.coolwarm, norm = None, 
                       frm = 'png', res = 300, sz = None, dim1 = True, dim2 = False, normalize = True,
                       return_handle = False, plt_handle = None, distance = 'euclidean', method = 'single', 
                       vmin = None, vmax = None, origin='lower'):

    if dim1 and mat.shape[0] > 0:
        print('Compute %s clustering in first dimension' % distance)
        if distance == 'nandist':
             dist = _nanDist(mat)
        else:
             dist = spsd.pdist(mat, distance)
             dist = spsd.squareform(dist)
     
        lnk1 = spch.linkage(dist, method=method, metric=distance)
#        dendro = spch.dendrogram(lnk1, p=500, no_plot=True, truncate_mode='mtica')
        dendro = spch.dendrogram(lnk1, p=100000, no_plot=True, truncate_mode=None)#{'mtica')
#        dendro = spch.dendrogram(lnk1, p=500, orientation='left', truncate_mode='mtica')
        lvs = dendro['leaves']
    else:
        lvs = np.arange(mat.shape[0])
        lnk1 = None
 
    if dim2 and mat.shape[0] > 1:
        print('Compute %s clustering in second dimension' % distance)
        if distance == 'nandist':
            dist = _nanDist(mat.T)
        else:
            dist = spsd.pdist(mat.T, distance)
            dist = spsd.squareform(dist)
 
        lnk2 = spch.linkage(dist, method=method, metric=distance)
#        dendro = spch.dendrogram(lnk2, p=500, no_plot=True, truncate_mode='mtica')
        dendro = spch.dendrogram(lnk2, p=100000, no_plot=True, truncate_mode=None)#'mtica')
        #dendro = spch.dendrogram(lnk2, p=100000, orientation='top', truncate_mode='mtica')
        lvs2 = dendro['leaves']
    else:
        lvs2 = np.arange(mat.shape[1])
        lnk2 = None

    cax = makeHeatmap(mat, fn=fn, tit=tit, xlab=xlab, ylab=ylab, cmap=cmap, norm=norm,
                      frm=frm, res=res, sz=sz, normalize=normalize, xidx=lvs, yidx=lvs2,
                      return_handle=True, plt_handle=plt_handle, vmin=vmin, vmax=vmax, origin=origin)
    if return_handle:
        return (cax, lvs, lvs2, lnk1, lnk2)
    else:
        return (lvs, lvs2, lnk1, lnk2)


def makeHeatmap(mat, fn = None, tit = None, xlab = None, ylab = None, cmap = cm.coolwarm, norm = None, 
               frm = 'png', res = 300, sz = None, normalize = True, xidx = None, yidx = None,
               return_handle = False, plt_handle = None, vmin = None, vmax = None, origin='lower'):

    if vmin is None:
        vmin = mat.min()
    if vmax is None:
        vmax = mat.max()

    if xidx is None:
        xidx = np.arange(mat.shape[0])
    if yidx is None:
        yidx = np.arange(mat.shape[1])
    
    if normalize and norm is None:
        norm = matplotlib.colors.Normalize(float(vmin)/2, float(vmax)/2) ##
  
    if not plt_handle:
        if sz is not None:
            fig = plt.figure(figsize = sz)
        else: 
            fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        ax = plt_handle

    if norm is not None:
        cax = ax.imshow(mat[xidx, :][:, yidx], interpolation='nearest', norm=norm, cmap=cmap, aspect='auto', origin=origin)
    else:
        cax = ax.imshow(mat[xidx, :][:, yidx], interpolation='nearest', cmap=cmap, aspect='auto', origin=origin)
    
    if tit is not None:
        ax.set_title(tit)
    if xlab is not None:
        ax.set_ylabel(xlab)
    if ylab is not None:
        ax.set_xlabel(ylab)
 
    if plt_handle is None:
        if not fn.endswith(frm):
             fn = "%s.%s" % (fn, frm)
        plt.savefig(fn, dpi=res, format=frm)
 
    if return_handle:
        return cax


def _nanDist(mat):
    """
    This is a helper function that inmplements a distance function ignoring
    NaN values in the input data. Not quite sure this works 100% ...
    """

    pdist = np.zeros((mat.shape[0], mat.shape[0]))
    for i in range(pdist.shape[0]):
        if i % 100 == 0:
            print('%i / %i' % (i, pdist.shape[0]))
        a = np.ones((pdist.shape[0] - i - 1, 1)) * mat[i, :] #, [pdist.shape[0] - i - 1, 1])
        b = mat[i + 1:pdist.shape[0], :]
        idx = (~np.isnan(a) * ~np.isnan(b))
        c = (a - b) * idx
        pdist[i, i+1:] = np.sum(c*c , axis = 1) / np.sum(idx, axis = 1)
    pdist += pdist.T
    pdist[np.isnan(pdist)] = np.nanmax(pdist)
    return pdist


def trackPlot(mat, fig=None, groups=None, ratios=None, labels=None, cmap=None, norm=None, is2D=False, xticks=False):
    """
    This function takes a matrix and generates a track figure with several panel according to a group structure that groups
    several rows/cols of the matrix into one panel. This can be done for rows only or for columns and rows. So if the input
    is a 10x10 matrix and we have a grouping of 2,4,3,1, then the final figure will have 4 panels, splitting the matrix into
    the respective groups. When option is2D is true, the same grouping is also applied to the columns. There is obviously room 
    for extension ...

    Input:
        mat - data matrix containing the values
        fig - figure object to place the panels into
        groups - grouping vector (is all ones per default, a single panel per row)
        ratios - the relative proportion each panel takes in the full plot (default to group values)
        labels - row labels for the matrix (needs to have as many entries as there are rows in the matrix)
        cmap - color map to apply to the single groups
        is2D - apply grouping to both columns and rows (rows only is default)
        xticks - set xticks

    Output:
        Returns a 2 tuple containing the figure object and an array with the axes objects corresponding to the single groups.
        fig - figure
        ax - axes
    """

    if fig is None:
        fig = plt.figure(figsize=(10, 10), dpi=200)
    if groups is None:
        groups = np.ones((mat.shape[0],), dtype='int')
    if ratios is None:
        ratios = groups
    if labels is not None:
        assert(labels.shape[0] == mat.shape[0])
    if cmap is None:
        cmap = np.array([plt.get_cmap('Blues')] * groups.shape[0], dtype='object')
    else:
        assert(cmap.shape[0] == groups.shape[0])
    if norm is None:
        norm = np.array([plt.Normalize(-1.0, 1.0)] * groups.shape[0], dtype='object')
    else:
        assert(norm.shape[0] == groups.shape[0])

    if is2D:
        gs = gridspec.GridSpec(groups.shape[0], groups.shape[0], height_ratios=ratios, hspace=0.05, width_ratios=ratios, wspace=0.05)
        last_col = 0
        axes = np.zeros((groups.shape[0], groups.shape[0]), dtype='object')
        for col in range(groups.shape[0]):
            last_row = 0
            for row in range(groups.shape[0]):
                axes[row, col] = fig.add_subplot(gs[row, col])
                axes[row, col].imshow(mat[last_row:last_row+groups[row], :][:, last_col:last_col+groups[col]], aspect='auto', origin='upper', interpolation='nearest', cmap=cmap[row], norm=norm[row])
                if xticks and row == 0:
                    axes[row, col].set_xticks(np.arange(groups[col]))
                    axes[row, col].xaxis.tick_top()
                    if labels is not None:
                        axes[row, col].set_xticklabels(labels[last_col:last_col+groups[col]], rotation=90)
                else:
                    axes[row, col].set_xticks([])
                if col == 0:
                    axes[row, col].set_yticks(np.arange(groups[row]))
                    if labels is not None:
                        axes[row, col].set_yticklabels(labels[last_row:last_row+groups[row]])
                else:
                    axes[row, col].set_yticks([])
                last_row += groups[row]
            last_col += groups[col]
    else:
        axes = np.zeros((groups.shape[0], ), dtype='object')
        gs = gridspec.GridSpec(groups.shape[0], 1, height_ratios=ratios, hspace=0.05)
        last_row = 0
        for row in range(groups.shape[0]):
            axes[row] = fig.add_subplot(gs[row, 0])
           # if density is not None and row in density:
           #     ax.fill_between(np.arange(mat.shape[1]), 
           # else:
            axes[row].imshow(mat[last_row:last_row+groups[row], :], aspect='auto', origin='lower', interpolation='nearest', cmap=cmap[row], norm=norm[row])
            axes[row].set_xticks([])
            axes[row].set_yticks(np.arange(groups[row]))
            if labels is not None:
                axes[row].set_yticklabels(labels[last_row:last_row+groups[row]])
            last_row += groups[row]

    return (fig, axes)

