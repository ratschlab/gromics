import sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy.random as npr
npr.seed(23)

import intervaltree as it
import scipy as sp
import scipy.stats as spst
import pdb
import h5py
import os
import re
import gzip
import pickle

from ratschlab_omicstools.viz.distribution import violin_plot

from spladder.classes import gene as cgene
from spladder.classes import splicegraph as csplicegraph
from spladder.classes import segmentgraph as csegmentgraph
sys.modules['modules.classes.gene'] = cgene
sys.modules['modules.classes.splicegraph'] = csplicegraph
sys.modules['modules.classes.segmentgraph'] = csegmentgraph

### command line parameters
if len(sys.argv) < 5:
    sys.stderr.write('Usage: %s <anno_junctions> <sample.count.hdf5> <sample.graph.pickle> <sample_junctions_w_outgroup.hdf5> <tcga_junctions> [<globsum>]\n' % sys.argv[0])
    sys.exit(1)

anno_junctions = sys.argv[1] # '/cluster/work/grlab/projects/TCGA/PanCancer/annotation/gencode.v19.annotation.hs37d5_chr.junctions.hdf5'
sample_count_file = sys.argv[2]
sample_event_file = sys.argv[3]
outgroup_junction_file = sys.argv[4]
tcga_neojunction_file  = sys.argv[5]

conf_count_glob = 20
if len(sys.argv) > 5:
    conf_count_glob = int(sys.argv[6])

### basic settings
CONF = 2
STRANDS = ['.', '+', '-']

conf_count = 3        # number of reads confirming a junction
gtex_thresh = 0.01    # fraction of samples in outgroup (GTEX) cohort for junction to be counted

counts_tcga = []    
counts_gtex = []
donor_support_tcga = []
log = False

# define output tags
conf_tag = '.globsum%i.conf%i' % (conf_count_glob, conf_count)
log_tag = ''
if log:
    log_tag = '.log10'
pickle_tag = '.G%s' % str(gtex_thresh) + conf_tag
plot_tag = log_tag + pickle_tag

outbase = re.sub(r'.hdf5$', '', sample_count_file)
sf_tc_pickle = outbase + '.sf.cpickle'
sf_gt_pickle = re.sub(r'.hdf5$', '', outgroup_junction_file) + '.sf.cpickle'

### preprocessing to harmonize genes used in the analysis
print('loading data for given sample')
IN_TC = h5py.File(sample_count_file, 'r')
gids_tcga = IN_TC['gene_ids_edges'][:, 0]
gnames_tcga = IN_TC['gene_names'][:].view(sp.chararray).decode('utf-8')
gid_names_tcga = sp.array([gnames_tcga[i] for i in gids_tcga], dtype='str') 
strains_tcga = IN_TC['strains'][:].view(sp.chararray).decode('utf-8')

lkidx = sp.arange(strains_tcga.shape[0])

### load sample event data
events, features = pickle.load(open(sample_event_file, 'rb'), encoding='latin1')

### load list of TCGA neojunctions
tcnj = sp.loadtxt(tcga_neojunction_file, dtype='str', delimiter='\t')
tcnj = sp.array([re.sub(':-:', ':m:', _) for _ in sp.unique(tcnj[:, 2])])
tcnj = sp.array([re.sub(r'-', ':', _) for _ in tcnj])
tcnj = sp.array([re.sub(':m:', ':-:', _) for _ in tcnj])
       
print('loading data for outgroup (GTEx)')
IN_GT = h5py.File(outgroup_junction_file, 'r')

if not os.path.exists(sf_gt_pickle):

    ### get gene intervals
    s_idx = sp.argsort(gids_tcga, kind='mergesort')
    _, f_idx = sp.unique(gids_tcga[s_idx], return_index=True)
    l_idx = sp.r_[f_idx[1:], gids_tcga.shape[0]]

    ### compute total edge count for outgroup samples
    print('Computing total edge count for outgroup (GTEx) samples')
    ### get counts
    genecounts_gtex = sp.zeros((f_idx.shape[0], IN_GT['edges_outgroup'].shape[1]), dtype='int')
    for i in range(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        genecounts_gtex[i, :] = sp.sum(IN_GT['edges_outgroup'][f_idx[i]:l_idx[i], :], axis=0)
    print('computing size factors for normalization')
    sf_gtex = spst.scoreatpercentile(genecounts_gtex, 75, axis=0)
    sf_gtex = sp.median(sf_gtex) / sf_gtex
    sf_gtex[sp.isnan(sf_gtex) | sp.isinf(sf_gtex)] = 1
    del genecounts_gtex

    ### store results in pickle
    pickle.dump(sf_gtex, open(sf_gt_pickle, 'wb'), -1)
else:
    print('loading size factors from %s' % sf_gt_pickle)
    sf_gtex = pickle.load(open(sf_gt_pickle, 'rb'))

if not os.path.exists(sf_tc_pickle):

    ### compute total edge count for given sample(s)
    print('Computing total edge count for given sample(s)')
    ### get gene intervals
    s_idx = sp.argsort(gids_tcga, kind='mergesort')
    _, f_idx = sp.unique(gids_tcga[s_idx], return_index=True)
    l_idx = sp.r_[f_idx[1:], gids_tcga.shape[0]]
    ### get counts
    genecounts_tcga = sp.zeros((f_idx.shape[0], lkidx.shape[0]), dtype='int')
    for i in range(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        genecounts_tcga[i, :] = sp.sum(IN_TC['edges'][s_idx[f_idx[i]:l_idx[i]], :][:, lkidx], axis=0)
    print('computing size factors for normalization')
    sf_tcga = spst.scoreatpercentile(genecounts_tcga, 75, axis=0)
    sf_tcga = sp.median(sf_tcga) / sf_tcga
    sf_tcga[sp.isnan(sf_tcga) | sp.isinf(sf_tcga)] = 1
    del genecounts_tcga

    ### store results in pickle
    pickle.dump(sf_tcga, open(sf_tc_pickle, 'wb'), -1)
else:
    print('loading size factors from %s' % sf_tc_pickle)
    sf_tcga = pickle.load(open(sf_tc_pickle, 'rb'))

compl_pickle = '%s%s.pickle' % (outbase, pickle_tag)
compl_pickle_junc = '%s%s_junc.hdf5' % (outbase, pickle_tag)
if not os.path.exists(compl_pickle) or not os.path.exists(compl_pickle_junc):

    ### remove edges that are confirmed with at least two reads in more than <gtex_thresh> (1%) of the outgroup (GTEx) samples
    print('Checking for edges that are common in provided outgroup (GTEx)')
    for i in range(0, IN_TC['edges'].shape[0], 1000):
        sys.stdout.write('.')
        if i > 0 and i % 10000 == 0:
            sys.stdout.write('%i/%i\n' % (i, IN_TC['edges'].shape[0]))
        sys.stdout.flush()
        tmp = sp.mean((IN_GT['edges_outgroup'][i:i+1000, :] * sf_gtex) >= 2, axis=1)
        if i == 0:
            k_idx = sp.where(tmp < gtex_thresh)[0]
        else:
            k_idx = sp.r_[k_idx, sp.where(tmp < gtex_thresh)[0] + i]
    print('Removed %i of %i edges that were common in outgroup (GTEx)' % (IN_TC['edges'].shape[0] - k_idx.shape[0], IN_TC['edges'].shape[0]))
    print('Retaining %i edges' % k_idx.shape[0])

    ### get coordinates of all junctions
    junction_coords = sp.c_[IN_GT['chrms'][:], IN_GT['strand'][:], IN_GT['pos'][:].astype('str')]
    junction_coords = sp.array([':'.join(x) for x in junction_coords])
    tcga_neojunction_match = sp.in1d(junction_coords, tcnj)

    ### remove edges that are not confirmed with at least conf_count_glob many reads in the current sample(s) and are not present in the TCGA neojunction list
    print('Checking for edges with insufficient support across samples')
    for i in range(0, IN_TC['edges'].shape[0], 1000):
        sys.stdout.write('.')
        if i > 0 and i % 10000 == 0:
            sys.stdout.write('%i/%i\n' % (i, IN_TC['edges'].shape[0]))
        sys.stdout.flush()
        tmp = (IN_TC['edges'][i:i+1000, :][:, lkidx] * sf_tcga).sum(axis=1)
        if i == 0:
            k_idx_ = sp.where((tmp > conf_count_glob) | tcga_neojunction_match[i:i+1000])[0]
        else:
            k_idx_ = sp.r_[k_idx_, sp.where((tmp > conf_count_glob) | tcga_neojunction_match[i:i+1000])[0] + i]
    kk_idx = sp.in1d(k_idx, k_idx_)
    print('Removed %i of %i edges that were not supported with at least %i reads in the sample cohort and are not present in the TCGA neojunction list' % (k_idx.shape[0] - sp.sum(kk_idx), k_idx.shape[0], conf_count_glob))
    print('Retaining %i edges' % sp.sum(kk_idx))
    k_idx = k_idx[kk_idx]

    ### remove edges that are annotated
    print('Checking for edges present in the annotation')
    INJA = h5py.File(anno_junctions, 'r')
    posj = sp.c_[INJA['chrms'][:], INJA['strand'][:], INJA['pos'][:].astype('str')]
    posj = sp.array([':'.join(x) for x in posj])
    INJA.close()

    kk_idx = sp.where(~sp.in1d(junction_coords[k_idx], posj))[0]
    print('Removed %i of %i edges that were found in the annotation' % (k_idx.shape[0] - kk_idx.shape[0], k_idx.shape[0]))
    print('Retaining %i edges' % kk_idx.shape[0])
    k_idx = k_idx[kk_idx]

    ### get gene intervals
    s_idx = sp.argsort(gids_tcga[k_idx], kind='mergesort')
    _, f_idx = sp.unique(gids_tcga[k_idx][s_idx], return_index=True)
    l_idx = sp.r_[f_idx[1:], k_idx.shape[0]]

    ### process GTEx samples
    print('\ncollecting counts for GTEx')

    ### get counts
    for i in range(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 500 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        counts_gtex.append(sp.sum((IN_GT['edges_outgroup'][sorted(k_idx[s_idx[f_idx[i]:l_idx[i]]]), :] * sf_gtex) >= conf_count, axis=0))

    ### process TCGA samples
    print('collecting counts for TCGA')
    junc_used = sp.zeros((k_idx.shape[0], lkidx.shape[0]), dtype='bool')

    ### get counts
    for i in range(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        tmp = IN_TC['edges'][sorted(k_idx[s_idx[f_idx[i]:l_idx[i]]]), :][:, lkidx] * sf_tcga
        counts_tcga.append(sp.sum(tmp >= conf_count, axis=0))
        donor_support_tcga.extend(sp.sum(tmp >= conf_count, axis=1))
        junc_used[s_idx[f_idx[i]:l_idx[i]], :] = (tmp >= conf_count)

    counts_tcga = sp.array(counts_tcga, dtype='int').T
    counts_gtex = sp.array(counts_gtex, dtype='int').T
    donor_support_tcga = sp.array(donor_support_tcga, dtype='int')

    pickle.dump((counts_tcga, counts_gtex, donor_support_tcga, k_idx), open(compl_pickle, 'wb'), -1)
    JU = h5py.File(compl_pickle_junc, 'w')
    JU.create_dataset(name='junc_used', data=junc_used, compression='gzip')
    JU.close()
else:
    print('loading counts from pickle')
    (counts_tcga, counts_gtex, donor_support_tcga, k_idx) = pickle.load(open(compl_pickle, 'rb'), encoding='latin1')
    JU = h5py.File(compl_pickle_junc, 'r')
    junc_used = JU['junc_used'][:]
    JU.close()

IN_TC.close()

### write junctions to file
jt_pos = IN_GT['pos'][:]
jt_chrm = IN_GT['chrms'][:].view(sp.chararray).decode('utf-8')
jt_strand = IN_GT['strand'][:].view(sp.chararray).decode('utf-8')
junc_out = gzip.open('%s%s_neojunctions.tsv.gz' % (outbase, pickle_tag), 'w')
for i in range(junc_used.shape[1]):
    strain = re.sub(r'.aligned', '', strains_tcga[i])
    iidx = sp.where(junc_used[:, i])[0]
    for ii in iidx:
        if len(gid_names_tcga.shape) > 1:
            gname = gid_names_tcga[k_idx[ii], 0]
        else:
            gname = gid_names_tcga[k_idx[ii]]
        txt = '\t'.join([strain, gname, '%s:%s:%i-%i' % (jt_chrm[k_idx[ii]], jt_strand[k_idx[ii]], jt_pos[k_idx[ii], 0], jt_pos[k_idx[ii], 1])]) + '\n'
        junc_out.write(txt.encode('utf-8'))
junc_out.close()
IN_GT.close()

