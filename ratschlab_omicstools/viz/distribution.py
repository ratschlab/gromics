import scipy as sp
from scipy.stats import gaussian_kde
import matplotlib
import matplotlib.pyplot as plt

def violin_plot(ax, data, pos, bp=False, fc='b'):
    '''
    create violin plots on an axis
    '''
 
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    i = -1
    for d,p in zip(data,pos):
        i += 1
        d_ = d[~sp.isnan(d)]
        if d_.shape[0] < 2:
            continue
        try:
            k = gaussian_kde(d_) #calculates the kernel density
        except:
            continue
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = sp.arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        _fc = fc[i] if isinstance(fc, list) else fc
        ax.fill_betweenx(x,p,v+p,facecolor=_fc,alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor=_fc,alpha=0.3)
    if bp:
        ax.boxplot(data,positions = pos)#,notch=1,positions=pos,vert=1)
    

def box_plot(ax, data, pos, fc='b', notch=False):
    '''
    create box plot using standardized framework
    '''

    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.75)
    i = -1
    for d,p in zip(data,pos):
        i += 1
        d_ = d[~sp.isnan(d)]
        if d_.shape[0] < 2:
            continue
        _fc = fc[i] if isinstance(fc, list) else fc
        medianprops=dict(linestyle='-', linewidth=3, color='k')
        capprops=dict(linestyle='-', color=_fc)
        boxprops=dict(linestyle='-', color=_fc)
        whiskerprops=dict(linestyle='--', color=_fc)
        ax.boxplot(d_, positions=[p], widths=[w], notch=notch, medianprops=medianprops, boxprops=capprops, capprops=capprops, whiskerprops=whiskerprops)


def qq_plot(pvals, ax=None, title=None, frm='png', dpi=150, fname=None, logscale=False):
    '''
    create a quantile quantile plot for the given p-values
    '''

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if title is not None:
        ax.set_title(title)
    ax.set_ylabel("Oberserved")
    ax.set_xlabel("Expected")

    exp = sp.linspace(0, 1, num=pvals.shape[0])
    if logscale and pvals.min() == 0:
        idx0 =  (pvals == 0)
        pvals[idx0] = pvals[~idx0].min() / 2

    if logscale:
        ax.plot(-sp.log10(exp), -sp.log10(sp.sort(pvals)), 'bo') 
    else:
        ax.plot(exp, sp.sort(pvals), 'bo') 
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ml = max(xlim[1], ylim[1])
    ax.plot([0, ml], [0, ml], 'r--')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if fname is not None:
        plt.savefig(fname, dpi=dpi, format=frm)

 
def dist_overview(data, ax=None, fname=None, format='pdf', log=False,
                  axis=0, sort=False):
    """
    Create a distribution overview plot

    data: measurements x samples
    ax: axis object (None)
    """

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if log:
        llq = sp.percentile(sp.log10(data + 1), 12.5, axis=axis)
        lq = sp.percentile(sp.log10(data + 1), 25, axis=axis)
        me = sp.percentile(sp.log10(data + 1), 50, axis=axis)
        uq = sp.percentile(sp.log10(data + 1), 75, axis=axis)
        uuq = sp.percentile(sp.log10(data + 1), 87.5, axis=axis)
    else:
        llq = sp.percentile(data, 12.5, axis=axis)
        lq = sp.percentile(data, 25, axis=axis)
        me = sp.percentile(data, 50, axis=axis)
        uq = sp.percentile(data, 75, axis=axis)
        uuq = sp.percentile(data, 87.5, axis=axis)

    if sort:
        s_idx = sp.argsort(me)
    else:
        s_idx = sp.arange(data.shape[axis])

    #ax.fill_between(sp.arange(data.shape[axis]), llq, uuq, facecolor='blue', edgecolor='blue', linestyle='--', alpha=0.2, interpolate=True)
    #ax.fill_between(sp.arange(data.shape[axis]), lq, uq, facecolor='blue', edgecolor='blue', linestyle='--', alpha=0.7, interpolate=True)
    ax.fill_between(sp.arange(data.shape[1-axis]), llq, uuq, facecolor='blue', edgecolor='none', alpha=0.2, interpolate=True)
    ax.fill_between(sp.arange(data.shape[1-axis]), lq, uq, facecolor='blue', edgecolor='none', alpha=0.7, interpolate=True)
    ax.plot(sp.arange(data.shape[1-axis]), me, 'r-')

    if fname is not None:
        plt.savefig(fname, format=format)

