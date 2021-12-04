"""
Useful plotting lib with various gwas standard plots
"""
import sys
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import scipy.stats as spst
import math

import pdb
import fdr


def _get_chrm_offsets(chr_lens):

    offsets = dict()
    s_keys = sorted(chr_lens)
    for i, k in enumerate(s_keys):
        if i == 0:
            offsets[k] = 0 
        else:
            offsets[k] = offsets[s_keys[i-1]] + chr_lens[s_keys[i-1]]
    return offsets

def _scalePosition(pos, chr_lens=None, N_chrm=0):
    """
    Returns a vector of relative position on the plot
    pos : Nx2 matrix of positions
    """

    offset    = 0
    if N_chrm == 0:
        N_chrm    = max(np.unique(pos[:,0]))#.shape[0]
    scaledPos = np.zeros(pos.shape[0])
    xTickPos  = np.zeros(N_chrm)

    chr_offsets = None
    if chr_lens is not None:
        chr_offsets = _get_chrm_offsets(chr_lens)

    for i,iChrm in enumerate(range(1,int(N_chrm)+1)):#np.unique(pos[:,0])):
#    for i, iChrm in enumerate(np.unique(pos[:,0])):
        #Get local stuff


        if chr_offsets is not None:
            offset = chr_offsets[iChrm]
        
        if (pos[:,0] == iChrm).sum() == 0: ### there are no points here
            if chr_offsets is not None:
                xTickPos[i] = offset + (0.5 * chr_lens[iChrm])
            continue
        locPos = pos[pos[:,0] == iChrm,:]
        
 
        #Scaling
        scaledPos[pos[:,0] == iChrm] = locPos[:,1] + offset
        if chr_offsets is None:
            xTickPos[i] = offset + (0.5 * max(locPos[:,1]))
            offset += max(locPos[:,1]) + 1
        else:
            xTickPos[i] = offset + (0.5 * chr_lens[iChrm])

    return scaledPos, xTickPos

def makeManhattanPlot(pv, pos, fn=None, qv=None, ax=None, gnPos=None, mutGene=None, tag = None):
    """
    Plots scales manhattan plot
    pos: Nx2 Matrix of chromsome and position
    pval: Nx1 Matrix of p-values
    gnPos: 3X1 chromosome, start, end. It is optional but if it is available it will plot an arrow at this position
    """

    ### check input consistency
    if np.sum(pv == 0) > 0:
        print('WARNING: Replacing p-Values with value 0 by 1!', file=sys.stderr)
        pv[pv==0] = 1
    if np.sum(np.isnan(pv)) > 0:
        print('WARNING: Replacing p-Values with value NaN by 1!', file=sys.stderr)
        pv[np.isnan(pv)] = 1

    ### compute FDR and number of contigs
    if qv is None:
        qv = fdr.qvalues(pv)
    else:
        assert(qv.shape[0] == pv.shape[0])
    Nchrm = max(np.unique(pos[:, 0]))#.shape[0]
    ### init variables
    scaledPos, xTickPos = _scalePosition(pos)

    ### Plotting
    if ax is None:
        plt.figure(figsize = [20,5])
        ax = plt.subplot(111)
    ax.set_xlim(min(scaledPos), max(scaledPos))
    ax.set_ylim(0, max(max(-np.log10(pv)), 10) + 1)
    ax.set_xlabel("Genomic Location")
    ax.set_ylabel("-log10(p-value)")
    if mutGene is not None:
        ax.set_title('MutGene: %s' % mutGene)
    ax.set_xticks(xTickPos)
    xtlab = ["Chr %i" %(x + 1) for x in range(Nchrm)]
    ax.set_xticklabels(xtlab, rotation=45)
    bnfThres = -np.log10(0.05 / pv.shape[0])
    colChr   = pos[:,0] % 2
 
    if np.sum(qv <= 0.05) != 0:
        #TODO: Make COlor dependent on Cancer Gene, TF or SF
        plt.plot(scaledPos[qv<=0.05], -np.log10(pv[qv<=0.05]), 'o', color='red')
    plt.scatter(scaledPos[qv>0.05], -np.log10(pv[qv>0.05]), c=colChr[qv>0.05], edgecolors='none')
    plt.plot([min(scaledPos), max(scaledPos)], [bnfThres, bnfThres], '--', linewidth=2, alpha=0.6)

    if gnPos is not None:
        iGnChr   = gnPos[0] == pos[:,0]
        iGnPosLB = gnPos[1] <= pos[:,1]
        iGnPosUB = gnPos[2] >= pos[:,1]
        if np.sum(iGnChr & iGnPosLB & iGnPosUB) != 0:
            gnStart  = scaledPos[iGnChr & iGnPosLB & iGnPosUB].min()
            gnEnd    = scaledPos[iGnChr & iGnPosLB & iGnPosUB].max()
            (y1, y2) = ax.get_ylim()
            plt.arrow(gnStart, y2 - 1, 0, -1 * math.ceil(y2 * 0.25), color='green', length_includes_head=True, head_length=0.4, width=math.ceil(0.05 * (max(scaledPos) - min(scaledPos))))

    if fn is not None:
        plt.savefig(fn, dpi = 100)


def makeManhattanPlot2D(pv, pos_gt, pos_pt, fn=None, ax=None, thresh=None, chrm_lens=None, title=None, color_chrms=False, marker_color='b', return_handle=False, label=None, return_pos = False, tag = None):
    """
    Plots scales manhattan plot
    pos_gt: Nx2 Matrix of chromsome and position for the genotype (p-values)
    pos_pt: Nx2 Matrix of chromsome and position for the phenotype (events/genes)
    pval: Nx1 Matrix of N p-values
    """

    ### check input consistency
    if np.sum(pv.ravel() == 0) > 0:
        print('WARNING: Replacing p-Values with value 0 by 1!', file=sys.stderr)
        pv[pv==0] = 1
    if np.sum(np.isnan(pv.ravel())) > 0:
        print('WARNING: Replacing p-Values with value NaN by 1!', file=sys.stderr)
        pv[np.isnan(pv)] = 1

    ### init variables
    Nchrm_gt = max(np.unique(pos_gt[:, 0]))#.shape[0]# np.unique(pos_gt[:,0]).shape[0]#
    Nchrm_pt = max(np.unique(pos_gt[:,0]))#.shape[0]#max(np.unique(pos_pt[:, 0]))#.shape[0]
    
    scaledPos_gt, TickPos_gt = _scalePosition(pos_gt, chrm_lens, N_chrm=22)
    scaledPos_pt, TickPos_pt = _scalePosition(pos_pt, chrm_lens, N_chrm=22)

    ### Plotting
    if ax is None:
        plt.figure(figsize = [10,10])
        ax = plt.subplot(111)
    if chrm_lens is None:
        ax.set_xlim(0, max(scaledPos_gt))
        ax.set_ylim(0, max(scaledPos_pt))
    else:
        ax.set_xlim(0, sum(chrm_lens.values()))
        ax.set_ylim(0, sum(chrm_lens.values()))
    ax.set_xlabel("Variant Position")
    ax.set_ylabel("Phenotype Position")
    ax.set_xticks(TickPos_gt)

    ax.set_xticklabels(["Chr %i" %(x + 1) for x in range(int(Nchrm_gt)+1)], rotation=45)
    ax.set_yticks(TickPos_pt)
    ax.set_yticklabels(["Chr %i" %(x + 1) for x in range(int(Nchrm_pt)+1)])
    if title is not None:
        ax.set_title(title)

    if color_chrms and chrm_lens is not None:
        chrm_offsets = _get_chrm_offsets(chrm_lens)
        s_keys = sorted(chrm_offsets)
        for k in range(2, len(s_keys), 2):
            ax.add_patch(Rectangle((chrm_offsets[s_keys[k-1]], 0), chrm_offsets[s_keys[k]] - chrm_offsets[s_keys[k-1]], ax.get_ylim()[1], color='grey', alpha=0.1))
            ax.add_patch(Rectangle((0, chrm_offsets[s_keys[k-1]]), ax.get_xlim()[1], chrm_offsets[s_keys[k]] - chrm_offsets[s_keys[k-1]], color='grey', alpha=0.1))
    if thresh is not None:
        k_idx = np.where(pv <= thresh)[0]
    else:
        k_idx = np.arange(pv.shape[0])

    if label is not None:
        if tag is None:
            h = ax.plot(scaledPos_gt[k_idx], scaledPos_pt[k_idx], 'o', color=marker_color, alpha='0.5', label=label) 
    else:
        if tag is None:
            h = ax.plot(scaledPos_gt[k_idx], scaledPos_pt[k_idx], 'o', color=marker_color, alpha=0.5)
        else:
            h = []
            for x in np.unique(tag):
                th, = ax.plot(scaledPos_gt[k_idx][tag[k_idx] == x], scaledPos_pt[k_idx][tag[k_idx] == x], 'o', alpha=0.2, label = x, markeredgecolor = 'none')
                h.append(th)

            ax.legend(loc = 'upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)

    if fn is not None:
        plt.savefig(fn, dpi = 100)

    if return_handle:
        return h
    if return_pos:
        return scaledPos_gt, scaledPos_pt
    

def simpleManhattan(l_pVal, l_pos, ax=None, fn=None, test=False, frm='pdf'):
    '''
    create a simple manhattan plot
    '''

    if fn is not None and not fn.endswith(frm):
        fn = '%s.%s' % (fn, frm)
    
    l_pos = np.array(l_pos, dtype = 'int')
    col = l_pos[:,0] / float(l_pos[:, 0].max())
 
    allChrm = np.unique(l_pos[:, 0])
    l_pVal = [-math.log10(float(x)) for x in l_pVal]
    bfr = -math.log10(0.05 / len(l_pos))
    
    if ax is None:
         fig = matplotlib.pyplot.figure()
         ax = fig.add_subplot(111)
    curmax = 0
    for rec in allChrm:
        l_pos[l_pos[:,0] == rec,1] += curmax
        curmax += l_pos[l_pos[:,0] == rec, 1].max()
    
    ax.scatter(l_pos[:,1], l_pVal, c = col)
    
    ax.set_ylabel("-log10(p)")
    ax.set_xlabel("Chromosome")
    ax.set_title("Manhattan Plot")

    if fn is not None:
        plt.savefig(fn, format=frm)


def evenSimplerManhattan(pv, pos, ax=None, fn=None):
    '''
    create a very very simple plot
    '''

    if ax is None:
        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)

    ax.plot(pos,-np.log10(pv),'.')[0]

    if fn is not None:
        plt.savefic(fn, format='pdf')

    if 'fig' in locals():
        plt.close(fig)
        


def qq_plot(pvals, ax=None, title=None, frm='png', res=150, fname=None, logscale=False):
    ''' 
    create a quantile quantile plot for the given p-values
    '''

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if title is not None:
        ax.set_title(title)

    exp = np.linspace(0, 1, num=pvals.shape[0] + 1)[1:]

    if logscale:
        ax.plot(-np.log10(exp), -np.log10(np.sort(pvals)), 'b.') 
        ax.set_ylabel("Observed (-log10)")
        ax.set_xlabel("Expected (-log10)")
    else:
        ax.plot(exp, np.sort(pvals)) 
        ax.set_ylabel("Observed")
        ax.set_xlabel("Expected")
    ax.plot([0, ax.get_xlim()[1]], [0, ax.get_xlim()[1]], 'g-')

    if fname is not None:
        plt.savefig(fname, dpi=res, format=frm)

    if 'fig' in locals():
        plt.close(fig)


