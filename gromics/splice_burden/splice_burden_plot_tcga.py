import sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np

color_dict = {'ACC': '#C1A72F', 'DLBC': '#3953A4', 'READ': '#DAF1FC', 'GBM': '#B2509E', 'THCA': '#F9ED32', 'BLCA': '#FAD2D9', 'UCEC': '#FBE3C7', 'PCPG': '#E8C51D', 'CESC': '#F6B667', 'UCS': '#F89420', 'THYM': '#CEAC8F', 'LIHC': '#CACCDB', 'CHOL': '#104A7F', 'HNSC': '#97D1A9', 'STAD': '#00AEEF', 'SKCM': '#BBD642', 'COAD': '#9EDDF9', 'UVM': '#009444', 'LUAD': '#D3C3E0', 'TGCT': '#BE1E2D', 'LUSC': '#A084BD', 'MESO': '#542C88', 'KIRC': '#F8AFB3', 'ESCA': '#007EB5', 'SARC': '#00A99D', 'KIRP': '#EA7075', 'LGG': '#D49DC7', 'PRAD': '#7E1918', 'PAAD': '#6E7BA2', 'BRCA': '#ED2891', 'OV': '#D97D25', 'KICH': '#ED1C24'}
log = False

def main():

    if len(sys.argv) < 6:
        sys.stderr.write('Usage: %s <neojunctions_sample> <neojunctions_tcga> <metadata_tcga> <sample_name> <outfname>\n' % sys.argv[0])
        sys.exit(1)

    neojunctions_file = sys.argv[1]
    neojunctions = np.loadtxt(neojunctions_file, dtype='str', delimiter='\t')

    tcga_neojunctions_file = sys.argv[2]
    tcga_neojunctions = np.loadtxt(tcga_neojunctions_file, dtype='str', delimiter='\t')

    tcga_metadata_file = sys.argv[3]
    tcga_metadata = np.loadtxt(tcga_metadata_file, dtype='str', delimiter='\t')
    tcga_metadata = tcga_metadata[1:, :]

    sample = sys.argv[4]
    outfname = sys.argv[5]

    ct_dict = dict([('%s.%s' % (_[2], _[4]), _[0]) for i,_ in enumerate(tcga_metadata)])
    tn_dict = dict([('%s.%s' % (_[2], _[4]), _[6] == 'False') for i,_ in enumerate(tcga_metadata)])

    ctypes = np.array([ct_dict[_] for _ in tcga_neojunctions[:, 0]])
    ctypes_u = np.unique(ctypes)
    tcga_is_tumor = np.array([tn_dict[_] for _ in tcga_neojunctions[:, 0]])

    tcga_neojunctions = tcga_neojunctions[:, 1].astype('int')

    ### plot sample within distribution
    fig = plt.figure(figsize=(18, 5), dpi=100)
    ax = fig.add_subplot(111)
    ### plot sums
    labels = []
    xticks = []
    cumx = 0
    buffsize = 200
    sort_means = []
    for ct in ctypes_u:
        tt_idx = np.where((ctypes == ct) & tcga_is_tumor)[0]
        sort_means.append(np.mean(tcga_neojunctions[tt_idx]))
    sort_means = np.array(sort_means)
    sidx = np.argsort(sort_means)

    for t in ctypes_u[sidx]:
        # tumor
        t_idx = np.where(ctypes == t)[0]
        if t_idx.size != 0:
            s_idx = np.argsort(tcga_neojunctions[t_idx])
            # tumor
            tt_idx = np.where(tcga_is_tumor[t_idx][s_idx])[0]
            if tt_idx.shape[0] != 0:
                ax.plot(np.arange(tt_idx.shape[0]) + cumx, tcga_neojunctions[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
                tn_idx = np.where(~tcga_is_tumor[t_idx][s_idx])[0]
                if tn_idx.shape[0] >= 5:
                    tn_median = np.median(tcga_neojunctions[t_idx[s_idx][tn_idx]])
                    ax.plot([cumx-25, t_idx.shape[0]+cumx+25], [tn_median, tn_median], ':r', linewidth=2.0)
            labels.append(t)
            xticks.append(cumx + int(t_idx.shape[0] / 2))
            cumx += t_idx.size + buffsize

    if log:
        ax.set_ylabel('Neojunctions (log10)')
    else:
        ax.set_ylabel('Neojunctions')
    #ax.set_xlabel('Samples'
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, rotation=90)
    ax.set_xlim([-1 * buffsize, cumx + buffsize])
    #ax.set_ylim([0, ymax])
    ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
    ax.xaxis.grid(False)

    ax.plot([0, ax.get_xlim()[1] - 25], [neojunctions.shape[0], neojunctions.shape[0]], 'b--', linewidth=2.5)
    ax.text(50, max(10, int(neojunctions.shape[0] * 1.01)), sample, color='b')

    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #axs.set_ticks_outer(ax)
    #axs.clean_axis(ax)
    ax.set_title('Splicing Complexity per cancer type')

    plt.tight_layout()
    outfname = outfname[:-4]
    plt.savefig(outfname + '.pdf', format='pdf')
    plt.savefig(outfname + '.png', format='png')
    plt.close(fig)

    return 0

if __name__ == "__main__":
    sys.exit(main())
