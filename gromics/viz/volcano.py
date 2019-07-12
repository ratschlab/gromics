"""
Volcano plotting lib
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')


def deseq_volcano(df, log2fc="log2FoldChange", qval="padj", labels=True, gene_name='external_gene_name', qval_threshold_line=0.01, l2fc_threshold_line=2, file_name="volcano_plot_rev.png"):
    """deseq csv volcano: import deseq csv file first as pandas dataframe
    df  dataframe   pandas dataframe
    log2fc  str column name that contains log2fc
    qval    str column name that contains padj
    labels  bool    label significant genes
    gene_name   str column name that contains gene name 
    qval_threshold_line  int    value to draw significant line
    l2fc_threshold_line  int    value to color points
    file_name   str where to save file
    """
    plt.close()
    fig, ax = plt.subplots()
    #color for significant points
    color_vec = []
    df['rev_fc'] = -1 * df[str(log2fc)]
    for idx in np.arange(df.shape[0]):
        if df.iloc[idx]['rev_fc'] > l2fc_threshold_line and df.iloc[idx][str(qval)] < 0.01:
            color_vec.append('r')
        elif df.iloc[idx]['rev_fc'] < -l2fc_threshold_line and df.iloc[idx][str(qval)] < 0.01:
            color_vec.append('b')
        else:
            color_vec.append('0.5')
    blue_cnt = color_vec.count('b')
    red_cnt = color_vec.count('r')
    ax.scatter(df['rev_fc'], -np.log10(np.array(df[str(qval)])), c=color_vec, alpha=0.5, edgecolors='none')
     
    #fix yaxis
    ax.set_ylim(-5, ax.get_ylim()[1])
    ax.set_ylabel('Significance (-log10 qvalue)')
    ax.set_xlabel('log2 Fold Change')
    ax.set_title('qvalue for Differential Expresison and log2 Fold Change\n%s downregulated (blue); %s upgregulated (red)\nqvalue < 0.01 and abs(log2FC) > 2' % (blue_cnt, red_cnt)) 
    plot_scale = ax.get_ylim()[1]
    scale = np.min([25, 0.2*plot_scale])
    scale = 25 
    if labels:
        df_sorted = df.sort_index(by=['padj']).iloc[0:25]
        for idx in np.arange(0, df_sorted.shape[0]):
            if idx % 4 == 0:
                plt.annotate(df_sorted.iloc[idx][str(gene_name)], \
                             xy = (df_sorted.iloc[idx]['rev_fc'], -np.log10(df_sorted.iloc[idx][str(qval)])), \
                             xytext = (-1*scale, scale), \
                             textcoords='offset points', \
                             ha='center', \
                             va='center', \
                             size='x-small', \
                             alpha=1, \
                             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            elif idx % 3 == 0:
                plt.annotate(df_sorted.iloc[idx][str(gene_name)], \
                             xy = (df_sorted.iloc[idx]['rev_fc'], -np.log10(df_sorted.iloc[idx][str(qval)])), \
                             xytext = (-1*scale, -1*scale), \
                             textcoords='offset points', \
                             ha='center', \
                             va='center', \
                             size='x-small', \
                             alpha=1, \
                             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            elif idx % 2 == 0:
                plt.annotate(df_sorted.iloc[idx][str(gene_name)], \
                             xy = (df_sorted.iloc[idx]['rev_fc'], -np.log10(df_sorted.iloc[idx][str(qval)])), \
                             xytext = (scale, scale), \
                             textcoords='offset points', \
                             ha='center', \
                             va='center', \
                             size='x-small', \
                             alpha=1, \
                             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            else:
                plt.annotate(df_sorted.iloc[idx][str(gene_name)], \
                             xy = (df_sorted.iloc[idx]['rev_fc'], -np.log10(df_sorted.iloc[idx][str(qval)])), \
                             xytext = (scale, -1*scale), \
                             textcoords='offset points', \
                             ha='center', \
                             va='center', \
                             size='x-small', \
                             alpha=1, \
                             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    plt.tight_layout()
    plt.savefig(file_name, dpi=200)
    plt.close()
